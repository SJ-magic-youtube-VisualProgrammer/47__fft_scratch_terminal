/************************************************************
source:
	C言語によるアルゴリズム辞典 p346
	
memo
	int fft(double x[], double y[], int IsReverse = FALSE);
		・標本点の数 "N"は2の階乗に限る
		・x[k]が実部、y[k]が虚部
		・結果はx[]、y[]に上書きする

参考URL
	■逆三角関数 : atan2
		https://ja.wikipedia.org/wiki/%E9%80%86%E4%B8%89%E8%A7%92%E9%96%A2%E6%95%B0
		
		■atan2
			https://ja.wikipedia.org/wiki/Atan2
		
	■窓関数使用時の補正！FFTの時に忘れがちな計算とは？
		https://watlab-blog.com/2019/04/20/window-correction/
	
	■窓関数を用いる理由
		https://www.logical-arts.jp/archives/124
		
	■窓関数の影響
		https://ecd-assist.com/wp/wp-content/uploads/2017/04/%E7%AA%93%E9%96%A2%E6%95%B0%E3%81%AE%E5%BD%B1%E9%9F%BF.pdf
************************************************************/

/************************************************************
************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/************************************************************
************************************************************/
static const int N = 64; // must be 2^x
static const double PI = 3.14159265398979;

double	sintbl[N + N/4];
int		bitrev[N];

float	fft_window[N];


/************************************************************
************************************************************/
int fft(double x[], double y[], int IsReverse = 0);
void make_bitrev(void);
void make_sintbl(void);


/************************************************************
************************************************************/

/******************************
******************************/
int main()
{
	/**************************
	初期化
	**************************/
	make_sintbl();
	make_bitrev();
	
	/**************************
	初期化 : 窓関数
	**************************/
	// 窓関数
	for(int i = 0; i < N; i++)	fft_window[i] = 0.5 - 0.5 * cos(2 * PI * i / N); // Hanning.
	
	// 補正係数
	float sum_window = 0;
	float acf; // Amplitude Correction Factor
	for(int i = 0; i < N; i++){
		sum_window += fft_window[i];
	}
	acf = sum_window / N;
	printf("> acf = %f\n", acf);
	fflush(stdout);
	
	/**************************
	x1[], y2[]に元となる波形を格納する
	**************************/
	int i;
	double x1[N], y1[N], x2[N], y2[N], x3[N], y3[N];

	for(i = 0; i < N; i++){
		// x1[i] = x2[i] = 6 * cos( (2*PI) * 3 * (i              ) / N );
		// x1[i] = x2[i] = 6 * cos( (2*PI) * 2 * (i + (double)N/8) / N );
		// x1[i] = x2[i] = 6 * cos( (2*PI) * 3 * (i + (double)N/8) / N );
		
		// x1[i] = x2[i] = 6 * cos( (2*PI) * 3 * (i) / N) + 4 * sin( (2*PI) * 9 * (i) / N);
		
		x1[i] = x2[i] = 6 * cos( (2*PI) * 3 * (i              ) / N );
		
		
		
		y1[i] = y2[i] = 0;
	}
	
	/**************************
	窓関数
	**************************/
	for(int i = 0; i < N; i++){
		x1[i] = x1[i] * fft_window[i];
		x2[i] = x2[i] * fft_window[i];
	}
	
	/**************************
	計算結果は、各時数のフーリエ係数 --- Excelで学ぶフーリエ変換で言う"Cn".
	x[]が実部、y[]が虚部

	Gain = 2 * sqrt(pow(x[], 2), pow(y[], 2));
	--- なぜ2倍が必要かはよくわからない.
	**************************/
	if(fft(x2, y2)) return 1;


	/**************************
	元に戻す.
	**************************/
	for(i = 0; i < N; i++){
		x3[i] = x2[i];
		y3[i] = y2[i];
	}
	
	if(fft(x3, y3, true)) return 1;


	/**************************
	Log.
	**************************/
	FILE* fp;
	fp = fopen("Log.csv", "w");
	if(fp == NULL){
		printf("File open Error\n");
		return 1;
	}
	for(i = 0; i < N; i++){
		fprintf(fp, "%d,%f,%f,,%d,%f,%f,,%d,%f,%f,,%d,%f\n", i, x1[i], y1[i], i, x2[i], y2[i], i, x3[i], y3[i], i, fft_window[i]);
	}
	fclose(fp);

	
	/**************************
	fin.
	**************************/
	printf("fin.\n");
	return 0;
}

/******************************
******************************/
int fft(double x[], double y[], int IsReverse)
{
	/*****************
		bit反転
	*****************/
	int i, j;
	for(i = 0; i < N; i++){
		j = bitrev[i];
		if(i < j){
			double t;
			t = x[i]; x[i] = x[j]; x[j] = t;
			t = y[i]; y[i] = y[j]; y[j] = t;
		}
	}

	/*****************
		変換
	*****************/
	int n4 = N / 4;
	int k, ik, h, d, k2;
	double s, c, dx, dy;
	for(k = 1; k < N; k = k2){
		h = 0;
		k2 = k + k;
		d = N / k2;

		for(j = 0; j < k; j++){
			c = sintbl[h + n4];
			if(IsReverse)	s = -sintbl[h];
			else			s = sintbl[h];

			for(i = j; i < N; i += k2){
				ik = i + k;
				dx = s * y[ik] + c * x[ik];
				dy = c * y[ik] - s * x[ik];

				x[ik] = x[i] - dx;
				x[i] += dx;

				y[ik] = y[i] - dy;
				y[i] += dy;
			}
			h += d;
		}
	}

	/*****************
	*****************/
	if(!IsReverse){
		for(i = 0; i < N; i++){
			x[i] /= N;
			y[i] /= N;
		}
	}

	return 0;
}

/******************************
******************************/
void make_bitrev(void)
{
	int i, j, k, n2;

	n2 = N / 2;
	i = j = 0;

	for(;;){
		bitrev[i] = j;
		if(++i >= N)	break;
		k = n2;
		while(k <= j)	{j -= k; k /= 2;}
		j += k;
	}
}

/******************************
******************************/
void make_sintbl(void)
{
	for(int i = 0; i < N + N/4; i++){
		sintbl[i] = sin(2 * PI * i / N);
	}
}
