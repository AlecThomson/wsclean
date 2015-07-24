#ifndef IMAGE_BUFFER_ALLOCATOR_H
#define IMAGE_BUFFER_ALLOCATOR_H

#include <complex>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <mutex>

//#define USE_DIRECT_ALLOCATOR

#ifndef USE_DIRECT_ALLOCATOR

template<typename NumType>
class ImageBufferAllocator
{
public:
	class Ptr
	{
	public:
		friend class ImageBufferAllocator<NumType>;
		Ptr()
			: _data(0) { }
		Ptr(NumType* data, ImageBufferAllocator<NumType>& allocator) 
			: _data(data), _allocator(&allocator) { }
		~Ptr()
		{
			reset();
		}
		NumType* data() const { return _data; }
		NumType& operator*() const { return *_data; }
		void reset()
		{
			if(_data != 0) _allocator->Free(_data);
			_data = 0;
		}
		void reset(NumType* data, ImageBufferAllocator<NumType>& allocator)
		{
			reset();
			_data = data;
			_allocator = &allocator;
		}
		NumType& operator[](size_t index) const { return _data[index]; }
	private:
		Ptr(const Ptr&) { }
		void operator=(const Ptr&) { }
		
		NumType* _data;
		ImageBufferAllocator<NumType>* _allocator;
	};
	
	ImageBufferAllocator() : _buffers(), _nReal(0), _nComplex(0), _nRealMax(0), _nComplexMax(0), _previousSize(0)
	{ }
	
	~ImageBufferAllocator()
	{
		std::lock_guard<std::mutex> guard(_mutex);
		size_t usedCount = 0;
		std::ostringstream str;
		for(typename std::vector<Buffer>::iterator i=_buffers.begin(); i!=_buffers.end(); ++i)
		{
			if(i->isFirstHalfUsed) {
				++usedCount;
				str << "Still used: buffer of " << i->size << '\n';
			}
			if(i->isSecondHalfUsed) {
				++usedCount;
				str << "Still used: buffer of " << i->size << '\n';
			}
			free(i->ptr);
		}
		if(usedCount != 0)
		{
			str << usedCount << " image buffer(s) were still in use when image buffer allocator was destroyed!";
			throw std::runtime_error(str.str());
		}
	}
	
	void ReportStatistics() const
	{
		std::lock_guard<std::mutex> guard(_mutex);
		double totalSize = 0.0;
		for(typename std::vector<Buffer>::const_iterator i=_buffers.begin(); i!=_buffers.end(); ++i)
		{
			totalSize += double(i->size) * double(sizeof(NumType)*2);
		}
		std::cout << "Image buf alloc stats:\n"
			"         max alloc'd images = " << _nRealMax << " real + " << _nComplexMax << " complex\n"
			"       max allocated chunks = " << _buffers.size() << "\n"
			"      current allocated mem = " << round(totalSize/1e8)/10.0 << " GB \n";
	}
	
	void Allocate(size_t size, Ptr& ptr)
	{
		ptr.reset(Allocate(size), *this);
	}
	
	NumType* Allocate(size_t size)
	{
		std::lock_guard<std::mutex> guard(_mutex);
		
		if(size != _previousSize)
		{
			FreeUnused();
			_previousSize = size;
		}
		
		++_nReal;
		if(_nReal > _nRealMax) _nRealMax = _nReal;
		for(typename std::vector<Buffer>::iterator i=_buffers.begin(); i!=_buffers.end(); ++i)
		{
			if(i->size == size)
			{
				if(!i->isFirstHalfUsed)
				{
					i->isFirstHalfUsed = true;
					return i->ptr;
				}
				if(!i->isSecondHalfUsed)
				{
					i->isSecondHalfUsed = true;
					return i->ptr + size;
				}
			}
		}
		Buffer* newBuffer = allocateNewBuffer(size);
		newBuffer->isFirstHalfUsed = true;
		return newBuffer->ptr;
	}
	
	std::complex<NumType>* AllocateComplex(size_t size)
	{
		std::lock_guard<std::mutex> guard(_mutex);
		
		if(size != _previousSize)
		{
			FreeUnused();
			_previousSize = size;
		}
		
		++_nComplex;
		if(_nComplex > _nComplexMax) _nComplexMax = _nComplex;
		for(typename std::vector<Buffer>::iterator i=_buffers.begin(); i!=_buffers.end(); ++i)
		{
			if(i->size == size)
			{
				if(!i->isFirstHalfUsed && !i->isSecondHalfUsed)
				{
					i->isFirstHalfUsed = true;
					i->isSecondHalfUsed = true;
					return reinterpret_cast<std::complex<NumType>*>(i->ptr);
				}
			}
		}
		Buffer* newBuffer = allocateNewBuffer(size);
		newBuffer->isFirstHalfUsed = true;
		newBuffer->isSecondHalfUsed = true;
		return reinterpret_cast<std::complex<NumType>*>(newBuffer->ptr);
	}
	
	void Free(NumType* buffer)
	{
		if(buffer != 0)
		{
			std::lock_guard<std::mutex> guard(_mutex);
			bool found = false;
			for(typename std::vector<Buffer>::iterator i=_buffers.begin(); i!=_buffers.end(); ++i)
			{
				if(i->ptr == buffer)
				{
					found = true;
					i->isFirstHalfUsed = false;
					break;
				}
				else if(i->ptr + i->size == buffer)
				{
					found = true;
					i->isSecondHalfUsed = false;
					break;
				}
			}
			if(!found)
				throw std::runtime_error("Invalid or double call to ImageBufferAllocator::Free(NumType*).");
			--_nReal;
		}
	}
	
	void Free(std::complex<NumType>* buffer)
	{
		if(buffer != 0)
		{
			std::lock_guard<std::mutex> guard(_mutex);
			bool found = false;
			for(typename std::vector<Buffer>::iterator i=_buffers.begin(); i!=_buffers.end(); ++i)
			{
				if(i->ptr == reinterpret_cast<NumType*>(buffer))
				{
					found = true;
					i->isFirstHalfUsed = false;
					i->isSecondHalfUsed = false;
					break;
				}
			}
			if(!found)
				throw std::runtime_error("Invalid or double call to ImageBufferAllocator::Free(std::complex<NumType>*).");
			--_nComplex;
		}
	}
	
	void FreeUnused()
	{
		size_t unusedCount = 0;
		typename std::vector<Buffer>::iterator i=_buffers.begin();
		while(i!=_buffers.end())
		{
			if(!i->isFirstHalfUsed && !i->isSecondHalfUsed)
			{
				free(i->ptr);
				_buffers.erase(i);
				i = _buffers.begin();
				++unusedCount;
			} else {
				++i;
			}
		}
		if(unusedCount != 0)
		{
			std::cout << "Freed " << unusedCount << " image buffer(s).\n";
		}
	}
	
private:
	struct Buffer
	{
		NumType* ptr;
		size_t size;
		bool isFirstHalfUsed, isSecondHalfUsed;
	};
	
	Buffer* allocateNewBuffer(size_t size)
	{
		_buffers.push_back(Buffer());
		Buffer* buffer = &_buffers.back();
		int errVal = posix_memalign(reinterpret_cast<void**>(&buffer->ptr), sizeof(NumType)*2, size * sizeof(NumType) * 2);
		if(errVal != 0)
		{
			switch(errVal)
			{
				case EINVAL:
					throw std::runtime_error("posix_memalign() failed: the alignment argument was not a power of two, or was not a multiple of sizeof(void *)");
				case ENOMEM:
					throw std::runtime_error("posix_memalign() failed: there was insufficient memory to fulfill the allocation request.");
				default:
					throw std::runtime_error("posix_memalign() failed but returned unknown error value");
			}
		}
		buffer->size = size;
		buffer->isFirstHalfUsed = false;
		buffer->isSecondHalfUsed = false;
		return buffer;
	}
	
	std::vector<Buffer> _buffers;
	size_t _nReal, _nComplex, _nRealMax, _nComplexMax, _previousSize;
	mutable std::mutex _mutex;
};

#else // USE_DIRECT_ALLOCATOR

template<typename NumType>
class ImageBufferAllocator
{
public:
	void ReportStatistics() const
	{
	}
	
	NumType* Allocate(size_t size)
	{
		return new NumType[size];
	}
	
	std::complex<NumType>* AllocateComplex(size_t size)
	{
		return new std::complex<NumType>[size];
	}
	
	void Free(NumType* buffer)
	{
		delete[] buffer;
	}
	
	void Free(std::complex<NumType>* buffer)
	{
		delete[] buffer;
	}
};

#endif // USE_DIRECT_ALLOCATOR

#endif
