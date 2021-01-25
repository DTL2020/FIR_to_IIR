#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum Weighting
{
    JINC,
    TRAPEZOIDAL
};

int iMul = 4;
int iTaps = 4;

//Weighting Weighting_type = TRAPEZOIDAL;//JINC;
Weighting Weighting_type = JINC;
float fPi = 3.14159265f;

float *g_pfKernel;
int iKernelSize;

float *input_x, *out_y;

float *coeff_b, *coeff_a;

int iN, iM, iL;

double *gauss_solver(double **a, double *y, int n) 
{
  double *x, max;
  int k, index;
  int i,j;
  const double eps = 0.000001;  // presision
  x = new double[n];
  k = 0;
  while (k < n) 
  {
    // Find line with max a[i][k]
    max = fabs(a[k][k]);
    index = k;
    for (i = k + 1; i < n; i++) 
    {
      if (fabs(a[i][k]) > max)
      {
        max = fabs(a[i][k]);
        index = i;
      }
    }
    // Reposition of lines
    if (max < eps) 
    {
/*      // not found any non-zero diagonal elements
      cout << "Can not got solution because of zero column ";
      cout << index << " marix A" << endl; */
      return 0;
    }
    for (j = 0; j < n; j++) 
    {
      double temp = a[k][j];
      a[k][j] = a[index][j];
      a[index][j] = temp;
    }
    double temp = y[k];
    y[k] = y[index];
    y[index] = temp;
    // Normalize equations
    for (i = k; i < n; i++) 
    {
      double temp = a[i][k];
      if (fabs(temp) < eps) continue; // For zero coefficient - skip
      for (j = 0; j < n; j++) 
        a[i][j] = a[i][j] / temp;
      y[i] = y[i] / temp;
      if (i == k)  continue; // Do not subtract equation from itself
      for (j = 0; j < n; j++)
        a[i][j] = a[i][j] - a[k][j];
      y[i] = y[i] - y[k];
    }
    k++;
  }
  // Backward setting
  for (k = n - 1; k >= 0; k--)
  {
    x[k] = y[k];
    for (i = 0; i < k; i++)
      y[i] = y[i] - a[i][k] * x[k];
  }
  return x;
}

float test_conv(void)
{
	// test with IIR convolution and measure error
	memset(input_x, 0, iL * 100 * sizeof(float));
	memset(out_y, 0, iL * 100 * sizeof(float));
	// 1-function input
    input_x[iL] = 1.0f;

	int i;

	for (i = 0; i <= iN; i++) 
	{
		printf(" coeff b%d = %f \r\n", i, coeff_b[i]);
	}
	for (i = 1; i <= iM; i++) 
	{
		printf(" coeff a%d = %f \r\n", i, coeff_a[i]);
	}

	float fMaxErr = 0;

	int n,k,l;

	for (n=0; n < iKernelSize/2; n++)
	{
		// b-coeff loop
		float fYn = 0.0f;

		for (k = 0; k <=iN; k++)
			fYn += coeff_b[k] * input_x[iL + n - k];

		// a-coeff loop
		for (l = 1; l <=iM; l++)
			fYn -= coeff_a[l] * out_y[iL + n - l];
		out_y[iL + n] = fYn;

		float fKv2 = g_pfKernel[(iKernelSize/2) * iKernelSize + (iKernelSize/2) + n];
		float fErr;

		if (fKv2 != 0)
		  fErr = (float)(fabs(fKv2 - fYn)/fabs(fKv2));
		else
			fErr = 0.0f;

		if (fErr > fMaxErr) fMaxErr = fErr;

		printf("Y %d = %f   fKv = %f   fErr = %f\r\n", n, fYn, fKv2, fErr);

	}

	return fMaxErr;
}

float test_conv_np(void)
{
	// test with IIR convolution and measure error
	memset(input_x, 0, iL * 100 * sizeof(float));
	memset(out_y, 0, iL * 100 * sizeof(float));
	// 1-function input
    input_x[iL] = 1.0f;

	float fMaxErr = 0;

	int n,k,l;

	for (n=0; n < iKernelSize/2; n++)
	{
		// b-coeff loop
		float fYn = 0.0f;

		for (k = 0; k <=iN; k++)
			fYn += coeff_b[k] * input_x[iL + n - k];

		// a-coeff loop
		for (l = 1; l <=iM; l++)
			fYn -= coeff_a[l] * out_y[iL + n - l];
		out_y[iL + n] = fYn;

		float fKv2 = g_pfKernel[(iKernelSize/2) * iKernelSize + (iKernelSize/2) + n];
		float fErr;

		if (fKv2 != 0)
		  fErr = (float)(fabs(fKv2 - fYn)/fabs(fKv2));
		else
			fErr = 0.0f;

		// limit max err at iN+iM to 0.05
		if ((n < iN+iM) && (fErr > 0.05f)) fErr = 100000.0f; 

		if (fErr > fMaxErr) fMaxErr = fErr;
	}

	return fMaxErr;
}


int main( void )
{
    int i, j;
	iKernelSize = iMul*iTaps*2;
	g_pfKernel = (float*)malloc(iKernelSize*iKernelSize*sizeof(float));
	memset(g_pfKernel, 0, iKernelSize*iKernelSize*sizeof(float));

    // our 2D just iMul-finely sampled kernel defined in output size


    /* Jinc weighted by Jinc - EWA Lanczos kernel */
    for (i = 0; i < iKernelSize; i++)
    {
        for (j = 0; j < iKernelSize; j++)
        {
            float fDist = sqrtf((float(iKernelSize / 2) - j) * (float(iKernelSize / 2) - j) + (float(iKernelSize / 2) - i) * (float(iKernelSize / 2) - i));

            // make kernel round in 2d
            if (fDist > iKernelSize / 2) continue;

            float fArg = (fPi * fDist / iMul);
//			fArg *= 1.1f;

            float fArg_w = fDist * 3.9f / (iKernelSize / 2);

            if (fArg != 0)
            {
                float fBess = 2.0f * (float)_j1(fArg) / fArg;
                float fW; // Jinc window
                if (Weighting_type == JINC)
                {
                     if (fArg_w != 0)
                     {
                         fW = 2.0f * (float)_j1(fArg_w) / fArg_w;
                     }
                     else
                     {
                         fW = 1.0f;
                     }
                }
                if (Weighting_type == TRAPEZOIDAL)
                {
                    fW = 1.0f;
                    if (fDist > (iKernelSize / 4))
                    {
                        fW = ((2 - (4 * fDist / iKernelSize))); // trapezoidal weighting
                    }
                }
                g_pfKernel[i * iKernelSize + j] = fBess * fW;
            }
            else
                g_pfKernel[i * iKernelSize + j] = 1.0f;


        } //j
    } //i


    // normalize to 1
    float fSum = 0.0f;
    for (i = 0; i < iKernelSize; i++)
    {
        for (j = 0; j < iKernelSize; j++)
        {
            fSum += g_pfKernel[i * iKernelSize + j];
        }
    }

    for (i = 0; i < iKernelSize; i++)
    {
        for (j = 0; j < iKernelSize; j++)
        {
            g_pfKernel[i * iKernelSize + j] /= fSum;
            g_pfKernel[i * iKernelSize + j] *= (iMul * iMul); // energy dissipated at iMul^2 output samples, so 1 norm * iMul^2
        }
    }


    for (i = iKernelSize/2; i < iKernelSize/2+1; i++)
    {
        for (j = 0; j < iKernelSize; j++)
        {
            float fKv = g_pfKernel[i * iKernelSize + j];
			printf ("i= %d , j = %d , fKv = %f\r\n", i,j, fKv); 

        }
    }

	iM = 3;
	iN = iM-1;
	iL = iM+iN+1;//10; // Num of equations in system

	input_x = new float[iL*100]; // centered to iL, x100 for future test
	out_y = new float[iL*100]; // centered to iL, x100 for future test

	memset(input_x, 0, iL * 100 * sizeof(float));
	memset(out_y, 0, iL * 100 * sizeof(float));

	// 1-function input
    input_x[iL] = 1.0f;

	for (i = 0; i < iL; i++)
	{
		out_y[i+iL] = g_pfKernel[(iKernelSize/2) * iKernelSize + (iKernelSize/2) + i];
	}

	double **a, *y, *solution;

	a = new double*[iL];
	y = new double[iL];

	int n,l,k;

	for (n = 0; n < iL; n++) // iL equations, y(n) each
	{
		a[n] = new double[iL];
        for (k = 0; k <= iN; k++)
		{
			a[n][k] = (double)input_x[iL + n - k]; 
		}

		for (l = 1; l <= iM; l++) 
		{
			a[n][iN+l] = (double) (-1.0 * out_y[iL + n - l]);
		}
	}

	for (i = 0; i < iL; i++) 
	{
		y[i] = g_pfKernel[(iKernelSize/2) * iKernelSize + (iKernelSize/2) + i]; // right part from center
		printf ("i= %d , y[i] = %f\r\n", i,y[i]); 
	}

	solution = gauss_solver(a, y, n);

//	coeff_b = new float[iN];
//	coeff_a = new float[iM+1];
	coeff_b = new float[iL];
	coeff_a = new float[iL];

	if (solution == 0)
	{
		printf ("solver return 0 error\r\n");
	}
	else
	{
		for (i = 0; i <= iN; i++) 
		{
			printf(" coeff b%d = %f \r\n", i, solution[i]);
			coeff_b[i] = (float)solution[i];
		}
		for (i = 1; i <= iM; i++) 
		{
			printf(" coeff a%d = %f \r\n", i, solution[i+iN]);
			coeff_a[i] = (float)solution[i+iN];
		}
	}

	// test with IIR convolution and measure error
	printf("initial solution\r\n");
//	test_conv();
	float fMinErr = test_conv();

	float c_a1_min_err = 0;
	float c_a2_min_err = 0;
	float c_a3_min_err = 0;
	float c_a4_min_err = 0;

	float c_b1_min_err = 0;
	float c_b2_min_err = 0;
	float c_b3_min_err = 0;

	float c_a2_init = coeff_a[2];
	float c_a3_init = coeff_a[3];
	float c_a4_init = coeff_a[4];

	float c_b1_init = coeff_b[1];
	float c_b2_init = coeff_b[2];
	float c_b3_init = coeff_b[3];

	coeff_a[1] -= coeff_a[1]*0.05f;

	int i1,i2,i3,i4,i5,i6,i7;

//	float fMulEP = 0.005f;
//	int iLps = 20;

	float fMulEP = 0.001f;
	int iLps = 100;

//	float fMulEP = 0.0001f;
//	int iLps = 1000;

	for (i1 = 0; i1 < iLps; i1++)
	{

		printf("outer loop %d of %d started\r\n", i1, iLps);

	    coeff_a[1] += coeff_a[1]*fMulEP;
		float fCurrErr = test_conv_np();
		if ( fCurrErr < fMinErr) 
		{
			fMinErr = fCurrErr;
			c_a1_min_err = coeff_a[1];
			c_a2_min_err = coeff_a[2];
			c_a3_min_err = coeff_a[3];
			c_a4_min_err = coeff_a[4];
			c_b1_min_err = coeff_b[1];
			c_b2_min_err = coeff_b[2];
			c_b3_min_err = coeff_b[3];
		}

		coeff_a[2] = c_a2_init - c_a2_init * 0.05f; 
		for (i2 = 0; i2 < iLps; i2++)
		{
			coeff_a[2] += coeff_a[2]*fMulEP;
			float fCurrErr = test_conv_np();
			if ( fCurrErr < fMinErr) 
			{
				fMinErr = fCurrErr;
				c_a1_min_err = coeff_a[1];
				c_a2_min_err = coeff_a[2];
				c_a3_min_err = coeff_a[3];
				c_a4_min_err = coeff_a[4];
				c_b1_min_err = coeff_b[1];
				c_b2_min_err = coeff_b[2];
				c_b3_min_err = coeff_b[3];
			}

		    coeff_a[3] = c_a3_init - c_a3_init * 0.05f; 
			for (i3 = 0; i3 < iLps; i3++)
			{
				coeff_a[3] += coeff_a[3]*fMulEP;
				float fCurrErr = test_conv_np();
				if ( fCurrErr < fMinErr) 
				{
					fMinErr = fCurrErr;
					c_a1_min_err = coeff_a[1];
					c_a2_min_err = coeff_a[2];
					c_a3_min_err = coeff_a[3];
					c_a4_min_err = coeff_a[4];
					c_b1_min_err = coeff_b[1];
					c_b2_min_err = coeff_b[2];
					c_b3_min_err = coeff_b[3];
				}

				coeff_b[1] = c_b1_init - c_b1_init * 0.05f; 
				for (i4 = 0; i4 < iLps; i4++)
				{
					coeff_b[1] += coeff_b[1]*fMulEP;
					float fCurrErr = test_conv_np();
					if ( fCurrErr < fMinErr) 
					{
						fMinErr = fCurrErr;
						c_a1_min_err = coeff_a[1];
						c_a2_min_err = coeff_a[2];
						c_a3_min_err = coeff_a[3];
						c_a4_min_err = coeff_a[4];
						c_b1_min_err = coeff_b[1];
						c_b2_min_err = coeff_b[2];
						c_b3_min_err = coeff_b[3];
					}

					coeff_b[2] = c_b2_init - c_b2_init * 0.05f; 
					for (i5 = 0; i5 < iLps; i5++)
					{
						coeff_b[2] += coeff_b[2]*fMulEP;
						float fCurrErr = test_conv_np();
						if ( fCurrErr < fMinErr) 
						{
							fMinErr = fCurrErr;
							c_a1_min_err = coeff_a[1];
							c_a2_min_err = coeff_a[2];
							c_a3_min_err = coeff_a[3];
							c_a4_min_err = coeff_a[4];
							c_b1_min_err = coeff_b[1];
							c_b2_min_err = coeff_b[2];
							c_b3_min_err = coeff_b[3];
						}
/*
						coeff_a[4] = c_a4_init - c_a4_init * 0.05f; 
						for (i6 = 0; i6 < iLps; i6++)
						{
							coeff_a[4] += coeff_a[4]*fMulEP;
							float fCurrErr = test_conv_np();
							if ( fCurrErr < fMinErr) 
							{
								fMinErr = fCurrErr;
								c_a1_min_err = coeff_a[1];
								c_a2_min_err = coeff_a[2];
								c_a3_min_err = coeff_a[3];
								c_a4_min_err = coeff_a[4];
								c_b1_min_err = coeff_b[1];
								c_b2_min_err = coeff_b[2];
								c_b3_min_err = coeff_b[3];
							}

							coeff_b[3] = c_b3_init - c_b3_init * 0.05f; 
							for (i7 = 0; i7 < iLps; i7++)
							{
								coeff_b[3] += coeff_b[3]*fMulEP;
								float fCurrErr = test_conv_np();
								if ( fCurrErr < fMinErr) 
								{
									fMinErr = fCurrErr;
									c_a1_min_err = coeff_a[1];
									c_a2_min_err = coeff_a[2];
									c_a3_min_err = coeff_a[3];
									c_a4_min_err = coeff_a[4];
									c_b1_min_err = coeff_b[1];
									c_b2_min_err = coeff_b[2];
									c_b3_min_err = coeff_b[3];
								}

							}

						}*/


					}

				}


			}

		}

	}

	printf("Min err %f b1 = %f b2 = %f a1 = %f a2 = %f a3 = %f\r\n", fMinErr, c_b1_min_err, c_b2_min_err, c_a1_min_err, c_a2_min_err, c_a3_min_err);

	coeff_a[1] = c_a1_min_err;
	coeff_a[2] = c_a2_min_err;
	coeff_a[3] = c_a3_min_err;
	coeff_a[4] = c_a4_min_err;

	coeff_b[1] = c_b1_min_err;
	coeff_b[2] = c_b2_min_err;
	coeff_b[3] = c_b3_min_err;

	printf("best solution\r\n");
	test_conv();


	free(g_pfKernel);


   return 0;
}