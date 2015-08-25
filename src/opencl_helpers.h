#ifndef OPEN_CL_HELPERS_H
#define OPEN_CL_HELPERS_H

#define SSS(string) #string

static int __cl_error = 0;

#define CL_CREATE_CONTEXT(context, device) \
context = clCreateContext(0, 1, &device, NULL, NULL, &__cl_error); \
if (!context) { \
printf("Error: Failed to create a compute context! %d\n", __cl_error); \
exit(1); \
}

#define CL_CREATE_IN_ORDER_COMMAND_QUEUE(queue, context, device) \
queue = clCreateCommandQueue(context, device, 0, &__cl_error); \
if (!queue) { \
printf("Error: Failed to create a command commands! %d\n", __cl_error); \
exit(1); \
}

#define CL_CREATE_PROGRAM_WITH_C_STR(program, context, c_str) \
program =  clCreateProgramWithSource(context, 1, &c_str, NULL, &__cl_error); \
if (!program) { \
printf("Error: Failed to create compute program! %d\n", __cl_error); \
exit(1); \
}

#define CL_BUILD_PROGRAM(program) \
if (clBuildProgram(program, 0, NULL, NULL, NULL, NULL) != CL_SUCCESS) {\
size_t len;\
char buffer[2048];\
printf("Error: Failed to build program executable!\n");\
clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);\
printf("%s\n", buffer);\
exit(1);\
}

#define CL_CREATE_KERNEL(kernel, program, kernel_name_str) \
kernel = clCreateKernel(program, kernel_name_str, &__cl_error); \
if (!kernel || __cl_error != CL_SUCCESS) { \
printf("Error: Failed to create compute kernel!\n"); \
exit(1);\
}

#define CL_CREATE_BUFFER(id, context, type, size) \
id = clCreateBuffer(context,  type,  size, NULL, NULL); \
if (!id) { \
printf("Error: Failed to allocate device memory!\n"); \
exit(1); \
}

#define CL_ENQUEUE_WRITE_BUFFER(queue, target_buffer, size, source_data) \
if (clEnqueueWriteBuffer(queue, target_buffer, CL_TRUE, 0, size, source_data, 0, NULL, NULL) != CL_SUCCESS) { \
printf("Error: Failed to write to source array!\n");\
exit(1);\
}

#define CL_SET_KERNEL_ARG(kernel, arg_position, size, source) \
__cl_error = clSetKernelArg(kernel, arg_position, size, & source); \
if (__cl_error != CL_SUCCESS) { \
printf("Error: Failed to set kernel arguments! %d\n", __cl_error);\
exit(1);\
}

#define CL_GET_LOCAL_WORK_GROUP_SIZE(kernel, device, size) \
__cl_error = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size), & size, NULL); \
if (__cl_error != CL_SUCCESS) { \
printf("Error: Failed to retrieve kernel work group info! %d\n", __cl_error); \
exit(1); \
}

#define CL_ENQUEUE_ND_RANGE_KERNEL(queue, kernel, dim, num_units, local_work_group_size) \
if (clEnqueueNDRangeKernel(queue, kernel, dim, NULL, & num_units, & local_work_group_size, 0, NULL, NULL)) { \
printf("Error: Failed to execute kernel!\n");\
exit(1); \
}

#define CL_ENQUEUE_ND_RANGE_KERNEL_AUTO_LOCAL(queue, kernel, dim, num_units) \
if (clEnqueueNDRangeKernel(queue, kernel, dim, NULL, & num_units, NULL, 0, NULL, NULL)) { \
printf("Error: Failed to execute kernel!\n");\
exit(1); \
}

#define CL_ENQUEUE_READ_BUFFER(queue, source_buffer, size, target_data) \
__cl_error = clEnqueueReadBuffer( queue, source_buffer, CL_TRUE, 0, size, target_data, 0, NULL, NULL ); \
if (__cl_error != CL_SUCCESS) { \
printf("Error: Failed to read output array! %d\n", __cl_error); \
exit(1); \
}

#define CL_GET_GPU_CONTEXT_AND_COMMAND_QUEUE(device, context, queue) \
if (clGetDeviceIDs(NULL, CL_DEVICE_TYPE_GPU, 1, &device, NULL) != CL_SUCCESS) { \
printf("Error: Failed to create a device group!\n"); exit(1); }\
CL_CREATE_CONTEXT(context, device);\
CL_CREATE_IN_ORDER_COMMAND_QUEUE(queue, context, device);

#define CL_GET_CPU_CONTEXT_AND_COMMAND_QUEUE(device, context, queue) \
if (clGetDeviceIDs(NULL, CL_DEVICE_TYPE_CPU, 1, &device, NULL) != CL_SUCCESS) { \
printf("Error: Failed to create a device group!\n"); exit(1); }\
CL_CREATE_CONTEXT(context, device);\
CL_CREATE_IN_ORDER_COMMAND_QUEUE(queue, context, device);

#define CL_CREATE_AND_BUILD_PROGRAM(program, context, source) \
CL_CREATE_PROGRAM_WITH_C_STR(program, context, source); \
CL_BUILD_PROGRAM(program);

struct __cl_releaser {
	void release(cl_mem m) {
		clReleaseMemObject(m);
	}
	void release(cl_command_queue q) {
		clReleaseCommandQueue(q);
	}
	void release(cl_program p) {
		clReleaseProgram(p);
	}
	void release(cl_kernel k) {
		clReleaseKernel(k);
	}
	void release(cl_context c) {
		clReleaseContext(c);
	}
	template <class T> __cl_releaser& operator,(T i) {
		release(i);
		return *this;
	}
};

#define CL_RELEASE(...) (__cl_releaser(),__VA_ARGS__);

#include <sstream>
#include <fstream>

using namespace std;

struct SourceReader {
	stringstream ss;
	string s;
	SourceReader() {}
	SourceReader(string filename) {
		ifstream t(filename);
		ss << t.rdbuf();
	}
	const char *getSource() {
		s = ss.str();
		return s.c_str();
	}
	operator const char *() {
		return getSource();
	}
	void insert(const string& s) {
		const string &temp = ss.str();
		ss.seekp(0);
		ss << s;
		ss << temp;
	}
	template <class T> void define(const string &name, T t) {
		const string &temp = ss.str();
		ss.seekp(0);
		ss << "#define " << name << " " << t << "\n";
		ss << temp;
	}
	void defineFloat(const string &name, float f) {
		const string &temp = ss.str();
		char output[32];
		sprintf(output, "%10f", f);
		ss.seekp(0);
		ss << "#define " << name << " " << output << "f\n";
		ss << temp;
	}
	
};


#endif