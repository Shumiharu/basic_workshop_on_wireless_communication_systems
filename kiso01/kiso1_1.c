#include <stdio.h>
#include <stdlib.h>		//乱数使用時のプログラム
#include <time.h>		//乱数の時間変化のプログラム
#include <math.h>		//数学的関数使用時のプログラム

#define PI 3.1415		//円周率の定義

int main(void)
{
	int i , j;   	
	int s;						//OS比
	double Ts;					//変調周期
	double f0;					//変調周波数
	double t=0.0;				//時間
	int a;						//波数
	int symble=8;   			//シンボル数
	double s_t;					//搬送波信号
	double IM,QM;				//I軸成分 , Q軸成分
	double P_base=0.0 ;			//ベースバンド信号の電力
	double P_carrier=0.0;		//搬送帯域信号の電力 
	
	///   OS比,周期，波数の入力    ///
	printf("OS比を入力してください : ");
	scanf("%d" , &s);
	
	printf("周期を入力してください : ");
	scanf("%lf" , &Ts);
	
	printf("Ts=%f\n", Ts*10000);
	printf("Ts=%f\n", Ts);
	
	
	printf("波数を入力してください : ");
	scanf("%d" , &a);
	
	//      周波数の計算       //
	f0=(double)a/Ts;

	printf("s=%d\n", s);
	printf("Ts=%f\n", Ts);
	printf("a=%d\n", a);
	printf("symble=%d\n", symble); 

	
	//      ランダム関数初期化      //
	srand((unsigned int)time(NULL));  	
	
	//		I軸、Q軸の信号サンプル		//
	int sI[s];		//I軸サンプル数
	int sQ[s];		//Q軸サンプル数
	
	//    1サンプル当たりの時間       //
	double T_sample;	//1サンプル当たりの時間	
	T_sample=Ts/(double)s;		
	
	//	   結果の表示のパラメータ　　　　//
	printf("   t          IM   　　 QM        s_t　\n");
	
	//  	ベースバンド信号I軸、Q軸及び搬送波信号の計算		//
	for(int i=0; i<symble; i++){	//1シンボルずつ計算
		sI[i]=rand()%2;			//送信ビット（I軸）
		sQ[i]=rand()%2;			//送信ビット（Q軸）
		
		for(int j=0; j<s; j++){			//1サンプルずつ計算
		if(sI[i]==0){
			if(sQ[i]==0)
			s_t=1/(sqrt(s))*cos(2*PI*f0*t+PI/4);		//送信ビット(0,0)
			else 
			s_t=1/(sqrt(s))*cos(2*PI*f0*t+3*PI/4);		//送信ビット(0,1)
		}	
			
		if(sI[i]==1){
			if(sQ[i]==1)
			s_t=1/(sqrt(s))*cos(2*PI*f0*t+5*PI/4);		//送信ビット(1,1)	
			else
			s_t=1/(sqrt(s))*cos(2*PI*f0*t+7*PI/4);		//送信ビット(1,0)
		}
		
	//     送信ビット0,1から１，-１に変換		//
		IM=1/(sqrt(2*s))*(-2*sI[i]+1);		//乱数が0ならプラス、1ならマイナスのベースバンド信号を表示
		QM=1/(sqrt(2*s))*(-2*sQ[i]+1);		//乱数が0ならプラス、1ならマイナスのベースバンド信号を表示
		
		printf("%f  ",t);			//時間tを表示
		printf("%f  " ,IM);     	//I軸の信号	                                                                                                                                                              );
		printf("%f  " ,QM);			//Q軸の信号
		printf("%f   \n",s_t);		//搬送波信号
	
		t=t+T_sample;				//次のシンボル時間に移行
		P_base=P_base+IM*IM+QM*QM;		//ベースバンド信号の電力の計算
		P_carrier=P_carrier+s_t*s_t;	//搬送帯域信号の電力の計算
		}
	}
	puts(" \n\n");
	printf("P_base=%f  ",P_base);				//ベースバンド信号の電力の表示
	printf("P_carrier=%f \n ",P_carrier);			//搬送波帯域信号の電力の表示
	return 0;
}
