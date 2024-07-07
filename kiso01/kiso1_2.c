#include <stdio.h>
#include <stdlib.h>		//乱数使用時のプログラム
#include <time.h>		//乱数の時間変化のプログラム
#include <math.h>		//数学的関数使用時のプログラム
#include <complex.h>	//複素数使用時のプログラム

#define PI 3.1415		//円周率の定義
#define Ts 0.001		//1シンボル時間
#define a 1				//波数
#define s 8				//オーバーサンプル比
#define symble 8		//シンボル数

int main(void)
{
	int i,j;  
	double f0;						//変調周波数
	double I_t[symble];				//複素ベースバンド信号の同相成分
	double Q_t[symble];				//複素ベースバンド信号の直交成分
	double complex u_t[symble];		//複素ベースバンド信号
	double s_t[symble];				//搬送帯域信号
	double t=0.0;					//時間
	double P_base=0.0;				//ベースバンド信号の電力
	double P_carrier=0.0;			//搬送帯域信号の電力
	
	
	//      周波数の計算       //
	f0=(double)a/Ts;
	
	//      ランダム関数初期化      //
	srand((unsigned int)time(NULL)); 
	
	//    1サンプル当たりの時間       //
	double T_sample;	//1サンプル当たりの時間	
	T_sample=Ts/(double)s;
	
	
	//	   結果の表示のパラメータ　　　　//
	printf("   t        I_t	　  Q_t      s_t　\n");

	for(i=0; i<symble; i++){
		I_t[i]=1/(sqrt(2*s))*(-2*(rand()%2)+1);	//電圧の大きさが1/√2s , 乱数が0なら１で1ならー1	
		Q_t[i]=1/(sqrt(2*s))*(-2*(rand()%2)+1);	//電圧の大きさが1/√2s , 乱数が0なら１で1ならー1
		
		u_t[i]=I_t[i]+Q_t[i]*I;		//u(t)=I(t)+jQ(t)
		
		for(j=0; j<s; j++){
			s_t[i]=creal(u_t[i])*cos(2*PI*f0*t)-cimag(u_t[i])*sin(2*PI*f0*t);
			
			printf("%f  ", t);
			printf("%f  ", creal(u_t[i]));
			printf("%f  ", cimag(u_t[i]));
			printf("%f  \n", s_t[i]);
			
			t=t+T_sample;
			P_base=P_base+I_t[i]*I_t[i]+Q_t[i]*Q_t[i];
			P_carrier=P_carrier+s_t[i]*s_t[i];
			
		}
	}
	printf("\n\n");
	printf("Pbase=%f W\n" , P_base);
	printf("Pcarrier=%f W\n" , P_carrier);
	return 0;
	
}

