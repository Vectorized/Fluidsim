#ifndef Vec_mm_vec_h
#define Vec_mm_vec_h

#define DEFAULT_SIZE 8
template <class T> class mini_vec {
	
private:
	
	size_t currentSize;
	size_t numElements;
	T *elements;
	
	inline size_t next_po2(size_t x)
	{
		int power = 2;
		while (x >>= 1) power <<= 1;
		return power;
	}
	
	inline void realloc(size_t size)
	{
		T * oldElements = elements;
		elements = new T[size];
		for (size_t i=0; i<numElements; ++i) elements[i] = oldElements[i];
		delete [] oldElements;
	}
	
	inline void resize()
	{
		if (numElements == currentSize) {
			currentSize <<= 1;
			realloc(currentSize);
		}
	}
	
public:
	
	inline void resize(size_t size)
	{
		if (size > numElements) {
			currentSize = next_po2(size);
			realloc(currentSize);
			numElements = size;
		}
	}
	
	mini_vec()
	{
		numElements = 0;
		currentSize = DEFAULT_SIZE;
		elements = new T[DEFAULT_SIZE];
	}
	
	inline bool empty()
	{
		return numElements == 0;
	}
	
	inline void push_back(const T &t)
	{
		elements[numElements++] = t;
		resize();
	}
	
	inline T pop_back()
	{
		T t = elements[--numElements];
		resize();
		return t;
	}
	
	inline bool contains(const T &t)
	{
		for (int i=0; i<numElements; ++i)
			if (elements[i] == t)
				return true;
		return false;
	}
	
	~mini_vec()
	{
		delete [] elements;
	}
	
	inline size_t size()
	{
		return numElements;
	}
	
	inline const T& operator [] (size_t i) const
	{
		return elements[i];
	}
	
	inline T& operator [] (size_t i)
	{
		return elements[i];
	}
	
	inline void clear()
	{
		numElements = 0;
	}
	
};
#undef DEFAULT_SIZE


#endif
