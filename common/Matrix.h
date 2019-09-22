#ifndef MATRIX_H
#define MATRIX_H

//#include <utility>
//#include <exception>
#include <memory>
#include <mex.h>
#include <matrix.h>

#define nullptr 0

using namespace std;

template<typename T> class Submatrix;
template<typename T> class WrapMatrix;

template<typename T>
class Matrix
{
public:
	friend class Submatrix<T>;
	friend class WrapMatrix<T>;

	Matrix()
		: m_iRowsNumber(0)
		, m_iColumnsNumber(0)
		, m_pElements(0)
		, m_bError(false)
	{
	}

	Matrix(size_t iRowsNumber, size_t iColumnsNumber, T* pElements = 0)
		: m_iRowsNumber(iRowsNumber)
		, m_iColumnsNumber(iColumnsNumber)
		, m_pElements(0)
		, m_bError(false)
		, m_iStrife(iColumnsNumber)
	{
		m_pElements = new T[m_iRowsNumber*m_iColumnsNumber];
		if(pElements)
		{
			std::memcpy(m_pElements, pElements, m_iRowsNumber*m_iColumnsNumber*sizeof(T));
		}
    }

	Matrix(const Matrix& matrix)
		: m_iRowsNumber(matrix.m_iRowsNumber)
		, m_iColumnsNumber(matrix.m_iColumnsNumber)
		, m_pElements(0)
		, m_bError(matrix.m_bError)
		, m_iStrife(matrix.m_iColumnsNumber)
	{
		if(m_pElements) delete[] m_pElements;
		m_pElements = new T[m_iRowsNumber*m_iColumnsNumber];
		if(matrix.m_pElements)
			if (matrix.m_iStrife == matrix.m_iColumnsNumber)
				memcpy(m_pElements, matrix.m_pElements, m_iRowsNumber*m_iColumnsNumber*sizeof(T));
			else
				for (size_t row=0; row<m_iRowsNumber; row++)
					memcpy(m_pElements + m_iStrife*row, matrix.getRowElements(row), m_iColumnsNumber*sizeof(T));
	}

	virtual ~Matrix()
	{
		if ( m_pElements ) delete[] m_pElements;
	}

	size_t getRowsNumber(void) const
	{
		return m_iRowsNumber;
	}

	size_t getColumnsNumber(void) const
	{
		return m_iColumnsNumber;
	}

	size_t getStrife(void) const
	{
		return m_iStrife;
	}

	const Matrix& operator=(const Matrix& matrix)
	{
		m_iRowsNumber = matrix.m_iRowsNumber;
		m_iColumnsNumber = matrix.m_iColumnsNumber;
		m_iStrife = matrix.m_iColumnsNumber;
		m_bError = matrix.m_bError;
		if(m_pElements)	delete[] m_pElements;
		m_pElements = new T[m_iRowsNumber*m_iColumnsNumber];
		if(matrix.m_pElements)
			if (matrix.m_iStrife == matrix.m_iColumnsNumber)
				memcpy(m_pElements, matrix.m_pElements, m_iRowsNumber*m_iColumnsNumber*sizeof(T));
			else
				for (size_t row=0; row<m_iRowsNumber; row++)
					memcpy(m_pElements + m_iStrife*row, matrix.getRowElements(row), m_iColumnsNumber*sizeof(T));
		return *this;
	}
	
	const T& operator[](std::pair<size_t, size_t> iIndex) const
	{
		return operator ()(iIndex.first, iIndex.second);
	}

	T& operator[](std::pair<size_t, size_t> iIndex)
	{
		return operator ()(iIndex.first, iIndex.second);
	}

	const T& operator()(std::pair<size_t, size_t> iIndex) const
    {
        return operator ()(iIndex.first, iIndex.second);
    }

    T& operator()(std::pair<size_t, size_t> iIndex)
    {
        return operator ()(iIndex.first, iIndex.second);
    }

	const T* operator[](int row) const
	{
		return m_pElements + row*m_iStrife;
	}

	T* operator[](int row)
	{
		return m_pElements + row*m_iStrife;
	}

	T& operator()(size_t iIndex) const
	{
		return m_pElements[iIndex];
	}

	T& operator()(size_t iRowIndex, size_t iColumnIndex)
	{
		return m_pElements[iRowIndex * m_iStrife + iColumnIndex];
	}

	const T& operator()(size_t iRowIndex, size_t iColumnIndex) const
	{
		return m_pElements[iRowIndex * m_iStrife + iColumnIndex];
	}

	T* getRowElements(size_t iRow)
	{
		return (m_pElements + iRow * m_iStrife);
	}

	const T* getRowElements(size_t iRow) const
	{
		return (m_pElements + iRow * m_iStrife);
	}

	T* getElements()
	{
		return m_pElements;
	}

	const T* getElements() const
	{
		return m_pElements;
	}

	Matrix operator+(Matrix& matrix) const
	{
		Matrix result(m_iRowsNumber, m_iColumnsNumber);
		if((m_iRowsNumber != matrix.m_iRowsNumber) || (m_iColumnsNumber != matrix.m_iColumnsNumber))
			result.m_bError = true;
		else
			for(size_t i = 0; i < m_iColumnsNumber; ++i)
				for(size_t j = 0; j < m_iColumnsNumber; ++j)
					result(i,j) = (*this)(i,j) + matrix(i,j);
		return result;
	}

	Matrix operator-(const Matrix& matrix) const
	{
		Matrix result(m_iRowsNumber, m_iColumnsNumber);
		if((m_iRowsNumber != matrix.m_iRowsNumber) || (m_iColumnsNumber != matrix.m_iColumnsNumber))
		{
			result.m_bError = true;
		}
		else
		{
			for(int i = 0; i < m_iRowsNumber*m_iColumnsNumber; ++i)
			{
				result(i) = (*this)(i) - matrix(i);
			}
		}
		return result;
	}

	Matrix operator*(const Matrix& matrix) const
	{
		Matrix result(m_iRowsNumber, matrix.m_iColumnsNumber);
		if(m_iColumnsNumber != matrix.m_iRowsNumber)
        {
			mexErrMsgTxt("ERROR: matrix dimensions must agree for multiplication");
        }
		else
		{
			for(size_t i = 0; i < result.m_iRowsNumber; ++i)
			{
				for(size_t k = 0; k < result.m_iColumnsNumber; ++k)
				{
					result(i, k) = (*this)(i, 0)*matrix(0,k);
					for(size_t j = 1; j < m_iColumnsNumber; ++j)
					{
						result(i, k) = result(i, k) + (*this)(i, j)*matrix(j, k); 
					}
				}
			}
		}
		return result;
	}

	Matrix transpose() const
	{
		Matrix result(m_iColumnsNumber, m_iRowsNumber);
		
		for(size_t i = 0; i < result.m_iRowsNumber; ++i)
		{
			for(size_t j = 0; j < result.m_iColumnsNumber; ++j)
			{
				result(i, j) = (*this)(j, i);
			}
		}

		return result;
	}

	bool isError() const
	{
		return m_bError;
	}

protected:
	size_t m_iRowsNumber;
	size_t m_iColumnsNumber;
	size_t m_iStrife;
	T*			 m_pElements;
	bool		 m_bError;
};

template<typename T>
class WrapMatrix : public Matrix<T>
{
public:
	WrapMatrix(size_t iRowsNumber, size_t iColumnsNumber, T* pElements)
		: Matrix<T>()
	{
		(*this).m_iRowsNumber = iRowsNumber;
		(*this).m_iColumnsNumber = iColumnsNumber;
		(*this).m_pElements = pElements;
		(*this).m_iStrife = iColumnsNumber;
	}

	virtual ~WrapMatrix()
	{
		(*this).m_pElements = nullptr;
	}
#ifdef MATLAB_MEX_FILE
	WrapMatrix(const mxArray* mat_array)
		: Matrix<T>()
    {
		(*this).m_iRowsNumber = mxGetN(mat_array);
		(*this).m_iColumnsNumber = mxGetM(mat_array);
		(*this).m_pElements = (T*)mxGetData(mat_array);
		(*this).m_iStrife = mxGetM(mat_array);
    }
#endif
};

template<typename T>
class Submatrix : public Matrix<T>
{
public:
	// FIXME: deal with const somehow
	Submatrix(const Matrix<T>& m, size_t first_row, size_t last_row, size_t first_column, size_t last_column)
		: Matrix<T>()
	{
		if (last_row > m.m_iRowsNumber)
			throw "Too many rows";
		if (last_column > m.m_iColumnsNumber)
			throw "Too many columns";
		(*this).m_iRowsNumber = last_row - first_row;
		(*this).m_iColumnsNumber = last_column - first_column;
		(*this).m_iStrife = m.m_iStrife;
		(*this).m_pElements = m.m_pElements + m.m_iStrife*first_row + first_column;
	}

	virtual ~Submatrix()
	{
		(*this).m_pElements = nullptr;
	}
};

#endif
