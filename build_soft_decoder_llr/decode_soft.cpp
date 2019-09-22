#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <errno.h>
#include <vector>
#include <list>

#include <mex.h>
#include <matrix.h>

#include "Matrix.h"
#include "common.h"
#include "LDPC.h"
#include "ObjectHandle.h"

#define ABS(A) (((A) >= 0) ? (A) : -(A))
using namespace std;

class LogLLR
{
public:
    double m_llr[Q];
    bool m_isnegative[Q];
public:
    LogLLR()
        : MAX_EXP(exp((double) 700))
        , MIN_EXP(exp((double)-700))
    {
        for(int i = 0; i < Q; ++i)
        {
            m_llr[i] = 0;
            m_isnegative[i] = false;
        }
    }

    LogLLR(const LogLLR& llr)
        : MAX_EXP(exp((double) 700))
        , MIN_EXP(exp((double)-700))
    {
        for(int i = 0; i < Q; ++i)
        {
            m_llr[i] = llr.m_llr[i];
            m_isnegative[i] = llr.m_isnegative[i];
        }
    }
    FieldElement getHard() const
    {
        double max = m_llr[0];
        int    index = 0;

        for (int i = 1; i < Q; ++i)
        {
            if (max <= m_llr[i])
            {
                max = m_llr[i];
                index = i;
            }
        }

        return FieldElement(index);
    }
    const LogLLR& operator=(const LogLLR& llr)
    {
        for(int i = 0; i < Q; ++i)
        {
            m_llr[i] = llr.m_llr[i];
            m_isnegative[i] = llr.m_isnegative[i];
        }
        return *this;
    }
    LogLLR operator*(const FieldElement& element) const
    {
        LogLLR result;
		FieldElement temp; 
			
        for(int i = 0; i < Q; ++i)
        {
			temp = FieldElement(i)*element;
            result.m_llr[temp.getElement()] = m_llr[i];
            result.m_isnegative[temp.getElement()] = m_isnegative[i];
        }

        return result;
    }
    const LogLLR& operator*=(const FieldElement& element)
    {
        LogLLR old = *this;
        FieldElement temp;

        for(int i = 0; i < Q; ++i)
        {
            temp = FieldElement(i)*element;
            m_llr[temp.getElement()] = old.m_llr[i];
            m_isnegative[temp.getElement()] = old.m_isnegative[i];
        }

        return *this;
    }
    LogLLR operator+(const LogLLR& llr) const
    {
        LogLLR result;

        for(int i = 0; i < Q; ++i)
        {
            result.m_llr[i] = m_llr[i] + llr.m_llr[i];
            result.m_isnegative[i] = (m_isnegative[i] != llr.m_isnegative[i]);
        }

        return result;
    }
    const LogLLR& operator+=(const LogLLR& llr)
    {
        for(int i = 0; i < Q; ++i)
        {
            m_llr[i] += llr.m_llr[i];
            m_isnegative[i] = (m_isnegative[i] != llr.m_isnegative[i]);
        }

        return *this;
    }
    LogLLR operator-(const LogLLR& llr) const
    {
        LogLLR result;

        for(int i = 0; i < Q; ++i)
        {
            result.m_llr[i] = m_llr[i] - llr.m_llr[i];
            result.m_isnegative[i] = (m_isnegative[i] != llr.m_isnegative[i]);
        }

        return result;
    }
    const LogLLR& operator-=(const LogLLR& llr)
    {
        for(int i = 0; i < Q; ++i)
        {
            m_llr[i] -= llr.m_llr[i];
            m_isnegative[i] = (m_isnegative[i] != llr.m_isnegative[i]);
        }

        return *this;
    }
    LogLLR operator-(double value) const
    {
        LogLLR result;

        for(int i = 0; i < Q; ++i)
        {
            result.m_llr[i] = m_llr[i] - value;
        }

        return result;
    }
    const LogLLR& operator-=(double value)
    {
        for(int i = 0; i < Q; ++i)
        {
            m_llr[i] -= value;
        }

        return *this;
    }
    void normalize()
    {
        double sum = 0;

        for (int i = 0; i < Q; ++i)
        {
            sum += safe_exp(m_llr[i]);
        }
        if (sum <= 0)
        {
            sum = Q;
        }
        sum = safe_log(sum);
        for (int i = 0; i < Q; ++i)
        {
            m_llr[i] -= sum;
        }
    }
    void sub_max()
    {
        int index = (*this).getHard().getElement();
        *this -= (*this).m_llr[index];
    }
    void sub_zero()
    {
        *this -= (*this).m_llr[0];
    }
    void fft() // Multi-dimensional Fourier Transform
    {
        for (int i = 0; i < Q; ++i)
        {
            m_llr[i] = safe_exp(m_llr[i]);
            if (m_isnegative[i])
            {
                m_llr[i] = -m_llr[i];
            }
        }
        ntt(m_llr);
        for (int i = 0; i < Q; ++i)
        {
            m_isnegative[i] = false;
            if (m_llr[i] < 0)
            {
                m_isnegative[i] = true;
            }
            m_llr[i] = safe_log(ABS(m_llr[i]));
        }
    }
    void print()
    {
        for (int i = 0; i < Q; ++i)
        {
            mexPrintf("%e ", m_llr[i]);
        }
        mexPrintf("\n");
    }
    void printExp()
    {
        for (int i = 0; i < Q; ++i)
        {
            mexPrintf("%e ", m_isnegative[i]? -safe_exp(m_llr[i]) : safe_exp(m_llr[i]));
        }
        mexPrintf("\n");
    }
private:
    const double MAX_EXP;
    const double MIN_EXP;
private:
    double safe_exp(double x)
    {
        if (x > 700)
        {
            return MAX_EXP;
        }
        if (x < -700)
        {
            return MIN_EXP;
        }
        return exp(x);
    }
    double safe_log(double x)
    {
        if (x < MIN_EXP)
        {
            return -700;
        }
        return log(x);
    }
    void ntt(double* p) // Multi-dimensional Fourier Transform
    {
        int factor = 1;
        for (int b = 0; b < LOG_Q; b++)
        {
            for (int rest = 0; rest < Q/2; rest++)
            {
                int restH = rest >> b;
                int restL = rest & (factor-1);
                int rest0 = (restH << (b+1)) + restL;
                int rest1 = rest0 + factor;
                double prest0 = p[rest0];
                p[rest0] += p[rest1];
                p[rest1] = prest0 - p[rest1];
            }
            factor += factor;
        }
    }
};

// Sum Product Decoder
bool SumProduct(LDPC& ldpc, vector<LogLLR>& in_llr, int max_iter, vector<LogLLR>& out_llr, vector<FieldElement>& y, int* number_of_iter)
{
    Matrix<LogLLR>  R_msgs(ldpc.m, ldpc.rmax);
    Matrix<LogLLR>  Q_msgs(ldpc.m, ldpc.rmax);
    Matrix<LogLLR>  Q_FFT(1, ldpc.rmax);

    for (int i = 0; i < ldpc.n; ++i)
    {
        out_llr[i] = in_llr[i];
        in_llr[i].sub_max();
	    y[i] = out_llr[i].getHard();
	  
        for (int j = 0; j < ldpc.col_weight[i]; ++j)
        {
            Q_msgs(ldpc.msgs_col(i,j)) = in_llr[i];
        }
    }

    for (int loop = 0; loop < max_iter; ++loop)
    {
        // Update R
        for (int j = 0; j < ldpc.m; j++)
        {
            LogLLR sum_FFT;
            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                Q_FFT(k) = Q_msgs(j, k)*ldpc.H(j,k);
                Q_FFT(k).fft();
                sum_FFT += Q_FFT(k);
            }
            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                R_msgs(j,k) = sum_FFT - Q_FFT(k);
                R_msgs(j,k).fft();
                R_msgs(j,k) = R_msgs(j,k)*(ldpc.H(j,k)^(-1));
            }
        }
        // Update Q
        for (int i = 0; i < ldpc.n; i++)
        {
            out_llr[i] = in_llr[i];
        
            for (int k = 0; k < ldpc.col_weight[i]; k++)
            {
                out_llr[i] += R_msgs(ldpc.msgs_col(i,k));
            }

            out_llr[i].sub_zero();
            y[i] = out_llr[i].getHard();
            
	        for (int k = 0; k < ldpc.col_weight[i]; k++)
		    {
		        Q_msgs(ldpc.msgs_col(i,k)) = out_llr[i] - R_msgs(ldpc.msgs_col(i,k));
		        Q_msgs(ldpc.msgs_col(i,k)).sub_max();
	        }
        }
    }

    if (number_of_iter)
    {
        *number_of_iter = max_iter;
    }
    return 1;
}

// Layered Sum Product Decoder
bool LayeredSumProduct(LDPC& ldpc, vector<LogLLR>& in_llr, int max_iter, vector<LogLLR>& out_llr, vector<FieldElement>& y, int* number_of_iter)
{
    Matrix<LogLLR>  R_msgs(ldpc.m, ldpc.rmax);
    Matrix<LogLLR>  Q_msg(1, ldpc.rmax);
    Matrix<LogLLR>  Q_FFT(1, ldpc.rmax);


    for (int i = 0; i < ldpc.n; ++i)
    {
        out_llr[i] = in_llr[i];
        y[i] = out_llr[i].getHard();
    }

    for (int loop = 0; loop < max_iter; ++loop)
    {
        // Update R
        for (int j = 0; j < ldpc.m; j++)
        {
            LogLLR sum_FFT;
            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                Q_msg(k) = out_llr[ldpc.row_col(j, k)] - R_msgs(j, k);
                Q_msg(k).sub_max();
                Q_FFT(k) = Q_msg(k)*ldpc.H(j,k);
                Q_FFT(k).fft();
                sum_FFT += Q_FFT(k);
            }
            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                R_msgs(j,k) = sum_FFT - Q_FFT(k);
                R_msgs(j,k).fft();
                R_msgs(j,k) = R_msgs(j,k)*(ldpc.H(j,k)^(-1));
                out_llr[ldpc.row_col(j, k)] = Q_msg(k) + R_msgs(j,k);
                out_llr[ldpc.row_col(j, k)].sub_zero();
                y[ldpc.row_col(j, k)] = out_llr[ldpc.row_col(j, k)].getHard();
            }
        }
    }

    if (number_of_iter)
    {
        *number_of_iter = max_iter;
    }
    return 1;
}

// Majority decoding
bool Majority(LDPC& ldpc, const vector<FieldElement>& x, int max_iter, vector<FieldElement>& y, int* number_of_iter)
{
    Matrix<FieldElement>  R_msgs(ldpc.m, ldpc.rmax);
    Matrix<FieldElement>  Q_msgs(ldpc.m, ldpc.rmax);
    
    for (int i = 0; i < ldpc.n; ++i)
    {
        y[i] = x[i];
        
        for (int j = 0; j < ldpc.col_weight[i]; ++j)
        {
            Q_msgs(ldpc.msgs_col(i,j)) = y[i];
        }
    }

    for (int loop = 0; loop < max_iter; ++loop)
    {
        // Update R
        for (int j = 0; j < ldpc.m; j++)
        {
            FieldElement sum;
            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                sum += Q_msgs(j, k)*ldpc.H(j,k);
            }
            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                R_msgs(j,k) = sum + Q_msgs(j, k)*ldpc.H(j,k); // Minus not defined
                R_msgs(j,k) = R_msgs(j,k)*(ldpc.H(j,k)^(-1));
            }
        }
        // Update Q
        for (int i = 0; i < ldpc.n; i++)
        {
            vector<int> result(Q, 0);
            for (int k = 0; k < ldpc.col_weight[i]; k++)
            {
                result[R_msgs(ldpc.msgs_col(i,k)).getElement()] += 1;
            }
            int max = 0;
            int index = -1;
            for (int k = 0; k < Q; k++)
            {
                if (result[k] >= max)
                {
                    max = result[k];
                    index = k;
                }
            }
            y[i] = index;
            for (int k = 0; k < ldpc.col_weight[i]; k++)
		    {
		        Q_msgs(ldpc.msgs_col(i,k)) = y[i];
	        }
        }
    }

    if (number_of_iter)
    {
        *number_of_iter = max_iter;
    }
    return 1;
}

// Majority Sequential decoding
bool MajoritySeq(LDPC& ldpc, const vector<FieldElement>& x, int max_iter, vector<FieldElement>& y, int* number_of_iter)
{
    Matrix<FieldElement>  R_msgs(ldpc.m, ldpc.rmax);
    Matrix<FieldElement>  Q_msgs(ldpc.m, ldpc.rmax);
    list<int> affected_rows;
    
    for (int i = 0; i < ldpc.n; ++i)
    {
        y[i] = x[i];
        
        for (int j = 0; j < ldpc.col_weight[i]; ++j)
        {
            Q_msgs(ldpc.msgs_col(i,j)) = y[i];
        }
    }

    for (int j = 0; j < ldpc.m; j++)
    {
        affected_rows.push_back(j);
    }

    for (int loop = 0; loop < max_iter; ++loop)
    {
        // Update Q
        for (int i = 0; i < ldpc.n; i++)
        {
            // Update R
            for (list<int>::const_iterator it = affected_rows.begin(); it != affected_rows.end(); ++it)
            {
                int j = *it;
                FieldElement sum;
                for (int k = 0; k < ldpc.row_weight[j]; k++)
                {
                    sum += Q_msgs(j, k)*ldpc.H(j,k);
                }
                for (int k = 0; k < ldpc.row_weight[j]; k++)
                {
                    R_msgs(j,k) = sum + Q_msgs(j, k)*ldpc.H(j,k); // Minus not defined
                    R_msgs(j,k) = R_msgs(j,k)*(ldpc.H(j,k)^(-1));
                }
            }
            affected_rows.clear();
            vector<int> result(Q, 0);
            for (int k = 0; k < ldpc.col_weight[i]; k++)
            {
                result[R_msgs(ldpc.msgs_col(i,k)).getElement()] += 1;
            }
            int max = 0;
            int index = -1;
            for (int k = 0; k < Q; k++)
            {
                if (result[k] >= max)
                {
                    max = result[k];
                    index = k;
                }
            }
            if (y[i].getElement() != index)
            {
                y[i] = FieldElement(index);
                for (int k = 0; k < ldpc.col_weight[i]; k++)
		        {
		            Q_msgs(ldpc.msgs_col(i,k)) = y[i];
	                affected_rows.push_back(ldpc.col_row(i,k));
                }
            }
        }
    }

    if (number_of_iter)
    {
        *number_of_iter = max_iter;
    }
    return 1;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
#if 0
    LogLLR a1;
    LogLLR a2;

    a1.m_llr[0] = log(1.0);
    a1.m_llr[1] = log(2.0);
    a1.m_llr[2] = log(3.0);
    a1.m_llr[3] = log(4.0);
    a1.m_llr[4] = log(5.0);
    a1.m_llr[5] = log(6.0);
    a1.m_llr[6] = log(7.0);
    a1.m_llr[7] = log(8.0);

    a1 -= a1.m_llr[0];
    a1.printExp();

    a2.m_llr[0] = log(2.0);
    a2.m_llr[1] = log(2.0);
    a2.m_llr[2] = log(3.0);
    a2.m_llr[3] = log(3.0);
    a2.m_llr[4] = log(4.0);
    a2.m_llr[5] = log(4.0);
    a2.m_llr[6] = log(5.0);
    a2.m_llr[7] = log(5.0);

    a2 -= a2.m_llr[0];
    a2.printExp();

    a1.fft();
    a1.printExp();

    a2.fft();
    a2.printExp();

    LogLLR a = a1 + a2;
    a.printExp();
    a.fft();
    a -= a.m_llr[0];

    a.printExp();

    return;
#endif
    int command = (int)*(mxGetPr(prhs[0]));
    switch(command)
    {
        case 0:
        {
            char buf[256];
            mxGetString(prhs[1], buf, 256); 
             
            LDPC* pLdpc = new LDPC();
            pLdpc->init(buf);
            
            if (nlhs > 0)
            {
	            plhs[0] = create_handle(pLdpc);
            }
            if (nlhs > 1)
            {
                plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
                *((double*)(mxGetPr(plhs[1]))) = pLdpc->n;
            }
            if (nlhs > 2)
            {	
                plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
                *((double*)(mxGetPr(plhs[2]))) = pLdpc->m;
            }
            break;   
        }
        case 1:
        case 2:
        {
            LDPC& ldpc = get_object<LDPC>(prhs[1]);

            vector<LogLLR> in_llrs(ldpc.n);
            vector<LogLLR> out_llrs(ldpc.n);
            vector<FieldElement> y(ldpc.n);

            double* llrs = mxGetPr(prhs[2]);
            int iMaxNumberOfIterations = (int)*(mxGetPr(prhs[3]));
            
            for(int i = 0; i < ldpc.n; ++i)
            {
                for(int j = 0; j < Q; ++j)
                {
                    in_llrs[i].m_llr[j] = llrs[i*Q + j];
                }
            }
         
            int number_of_iter = 0;
         
            bool result = true;
            if (command == 1)
            {
                result = SumProduct(ldpc, in_llrs, iMaxNumberOfIterations, out_llrs, y, &number_of_iter);
            }
            else if (command == 2)
            {
                result = LayeredSumProduct(ldpc, in_llrs, iMaxNumberOfIterations, out_llrs, y, &number_of_iter);
            }
                     
            /* denial flag */
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[0]))) = result;
         
            /* number of iterations */
		    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
		    *((double*)(mxGetPr(plhs[1]))) = number_of_iter;
		 
		    double* data = NULL;
		 
		     /* hard desicion */
		     plhs[2] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
		     data = mxGetPr(plhs[2]);
		     for(int i = 0; i < ldpc.n; ++i)
		     {
			     data[i] = y[i].getElement();
		     }
    		 
		     /* a posteriori llr */
		     plhs[3] = mxCreateDoubleMatrix(Q, ldpc.n, mxREAL);
		     data = mxGetPr(plhs[3]);
		     for(int i = 0; i < ldpc.n; ++i)
		     {
			     for(int j = 0; j < Q; ++j)
			     {
				     data[i*Q + j] = out_llrs[i].m_llr[j];
			     }
		     }
		     break;
        }
        case 3:
        case 4:
        {
            LDPC& ldpc = get_object<LDPC>(prhs[1]);

            vector<FieldElement> x(ldpc.n);
            vector<FieldElement> y(ldpc.n);

            double* input = mxGetPr(prhs[2]);
            int iMaxNumberOfIterations = (int)*(mxGetPr(prhs[3]));
            
            for(int i = 0; i < ldpc.n; ++i)
            {
                x[i] = FieldElement((int)input[i]);
            }
         
            int number_of_iter = 0;
         
            bool result = true;
            if (command == 3)
            {
                result = Majority(ldpc, x, iMaxNumberOfIterations, y, &number_of_iter);
            }
            else if (command == 4)
            {
                result = MajoritySeq(ldpc, x, iMaxNumberOfIterations, y, &number_of_iter);
            }
                     
            /* denial flag */
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[0]))) = result;
         
            /* number of iterations */
		    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
		    *((double*)(mxGetPr(plhs[1]))) = number_of_iter;
		 
		    double* data = NULL;
		 
		     /* hard desicion */
		     plhs[2] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
		     data = mxGetPr(plhs[2]);
		     for(int i = 0; i < ldpc.n; ++i)
		     {
			     data[i] = y[i].getElement();
		     }
    		 
		     break;
        }
        default:
        break;
    }
}

