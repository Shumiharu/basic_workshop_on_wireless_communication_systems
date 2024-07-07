#include <stdio.h>
#include <stdlib.h>		//乱数使用時のプログラム
#include <time.h>		//time関数(乱数の種を設定するため)使用時のプログラム
#include <math.h>		//数学的関数使用時のプログラム
#include <complex.h>	//複素数使用時のプログラム

#define Ts 0.001				//1シンボル時間
#define a 1						//波数
#define s 4						//オーバーサンプル比
#define symble 1024				//シンボル数
#define f0 (a/Ts)				//変調周波数
#define T_sample (Ts/s) 		//1サンプル当たりの時間
#define N (s*symble)			//全サンプル数
#define delta_f (1/(Ts*symble))	//Δf=1/T Δf:スペクトルの幅  T(=Ts*symble):時間長
#define loop 100				//試行回数
//確率密度関数用↓↓↓↓
#define range 1				//試行範囲
#define w 0.05					//範囲の幅


//BER特性用↓↓↓↓
#define c 2 					//1シンボル当たりのビット数 QPSKの場合は2
#define BnTs 1    				//1シンボル当たりb点を用いるとき，BnTs=b

/////////////////////////　　　　FFT　　　　//////////////////////////////
void FFT(double complex Time[] , double complex Freq[]){

	int i,m,k,j;
	double arg;
	double complex X[N],Y[N];
	
	for(i=0;i<N; i++){
		Freq[i]=0;
		X[i]=Time[i];
	}
	for(m=2; m<=N; m*=2){		//バタフライの段数
		j=m/2;		//何個とびかを表す
		
		for(i=0; i<N/m; i++){		//
			for(k=0; k<j; k++){			//どこの数字と計算するか，何個とびか
				arg=2*M_PI*(k-N/2)/m;
				Y[i*m+k]=X[i*j+k]+cexp(-arg*I)*X[i*j+(N/2)+k];
				Y[i*m+j+k]=X[i*j+k]-cexp(-arg*I)*X[i*j+(N/2)+k];
			}
		}
		for(i=0; i<N; i++){
			X[i]=Y[i];
		}
		for(i=0; i<N; i++)
			Y[i]=0.0+0.0*I;
	}
	
	for(i=0;i<N;i++)
		Freq[i]=X[i]/sqrt((double)N);
}

/////////////   複素ベースバンド信号(と搬送帯域信号)の生成　　　//////////////////////
///　　搬送帯域信号も欲しいなら☆の部分のコメント化を解除する   ////
// void baseband_carrier(const double time[] , double complex base[] , double complex carrier[]){
// 		double tmp_I;	//同相成分Iを一時的に格納
// 		double tmp_Q;	//直交成分Qを一時的に格納
		
// 		for(int i=0; i<symble; i++){
// 			tmp_I=1/(sqrt(2*s))*(-2*(rand()%2)+1);	//電圧の大きさが1/√2s , 乱数が0なら１で1ならー1	
// 			tmp_Q=1/(sqrt(2*s))*(-2*(rand()%2)+1);	//電圧の大きさが1/√2s , 乱数が0なら１で1ならー1
	
// 			for(int j=0; j<s; j++){
// 				base[s*i+j]=tmp_I+tmp_Q*I;				//複素ベースバンド信号u(t)=I(t)+jQ(t)
// 				//☆ carrier[s*i+j]=creal(base[s*i+j])*cos(2*M_PI*f0*time[s*i+j])-cimag(base[s*i+j])*sin(2*M_PI*f0*time[s*i+j]);	//搬送帯域信号s(t)
			
// 				// printf("%f  ", time[s*i+j]);
// 				// printf("%f  ", creal(base[s*i+j]));
// 				// printf("%f  ", cimag(base[s*i+j]));
// 				// printf("%f  \n", carrier[s*i+j]);
// 			}
// 		}	
// }



//　　　　Box-Muller法によるガウス雑音生成　　　　　　//
void boxmuller(double complex noise[], double mean , double variance){
  int i;
	double u1,u2;
	double nreal[N],nimag[N];
	
	for(i=0; i<N; i++){
	////////////////　　　　[0,1]の一様乱数生成　　　　　/////////////

		u1=((double)rand()+1.0)/((double)RAND_MAX+1.0);		//rand()=0～32767
		u2=((double)rand()+1.0)/((double)RAND_MAX+1.0);
		
		///   標準偏差σを1/√2にすれば各軸の分散σ^2=1/2になる．
		nreal[i]=sqrt(variance)*sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2)+mean;		//分散を1/√2にして合成後の電力を1Wにする．
		nimag[i]=sqrt(variance)*sqrt(-2.0*log(u1))*sin(2.0*M_PI*u2)+mean;
		noise[i]=nreal[i]+nimag[i]*I;
	}
}

////////////    BER_QPSKを求める関数　　　　　/////////////　　　		　　　　　　
void BER_QPSK(double time[], double variance , int rho){
	double ci,cq;						//BERカウント用
	double cs;							//SERカウント用
	double complex base_QPSK[N];		//QPSKの複素ベースバンド信号
	double complex carrier_QPSK[N];		//QPSKの搬送帯域信号（形式上)
	double complex awgn_QPSK[N];		//複素ガウス雑音
	double complex receiver[N];			//QPSKの受信信号r(t)
	long attempt=0;						//試行回数
	double ber,ser;						//誤り率
	
	ci=0.0;	cq=0.0; cs=0.0;
	//////※ρ=9とか10以上の時はdo{ }whileでciとcqの個数でループさせる////////////
	// do{
		/*
	for(int m=0; m<loop; m++){
		baseband_carrier(time,base_QPSK,carrier_QPSK);	//loop毎にベースバンド信号生成
		boxmuller(awgn_QPSK,0.0,variance);				//loop毎にガウス雑音生成
		
		for(int i=0; i<N; i++)
			receiver[i]=base_QPSK[i]+awgn_QPSK[i];		//受信信号r(t)=u(t)+n(t)
			
		for(int j=0; j<N; j++){
			////////////  BER  ////////////
			if(creal(base_QPSK[j])>=0 && creal(receiver[j])<0) ci++;
			if(creal(base_QPSK[j])<0 && creal(receiver[j])>=0) ci++;
			if(cimag(base_QPSK[j])>=0 && cimag(receiver[j])<0) cq++;	
			if(cimag(base_QPSK[j])<0 && cimag(receiver[j])>=0) cq++; 
			
/* 			////////////  SER   ///////////
			if(creal(base_QPSK[j])>=0 && cimag(base_QPSK[j])>=0){	//信号点が第一象限の時
				if(creal(receiver[j])<0 || cimag(receiver[j])<0) cs++;
			}
			if(creal(base_QPSK[j])<0 && cimag(base_QPSK[j])>=0){	//信号点が第二象限の時
				if(creal(receiver[j])>=0 || cimag(receiver[j])<0) cs++;
			}
			if(creal(base_QPSK[j])<0 && cimag(base_QPSK[j])<0){		//信号点が第三象限の時
				if(creal(receiver[j])>=0 || cimag(receiver[j])>=0) cs++;
			}
			if(creal(base_QPSK[j])>=0 && cimag(base_QPSK[j])<0){	//信号点が第四象限の時
				if(creal(receiver[j])<0 || cimag(receiver[j])>=0) cs++;
			} 
		}
		attempt++;
	}*/
	// }while((ci+cq)<=50000.0);	//カウントの合計が---になったら終了
	///////  送信信号と受信信号の表示　　/////////
/* 	if(rho==10){
		printf("　　　　送信信号u(t)と受信信号r(t)の表示　　　\n");
		printf("　時間t　　       u(t)のI軸         u(t)のQ軸             r(t)のI軸              r(t)のQ軸 \n"); 
		for(int i=0; i<N; i++){
			printf("%f          " , time[i]);
			printf("%f          " , creal(base_QPSK[i]));
			printf("%f          " , cimag(base_QPSK[i]));
			printf("%.10f          " , creal(receiver[i]));		
			printf("%.10f\n" , cimag(receiver[i]));
		}	
	} */
	
	///////  BERの計算と表示   //////////
	// printf("ci=%f , cq=%f \n" , ci,cq);
	ber=(ci+cq)/(2.0*(double)attempt*N);
	printf(" %d          %.10f \n",rho,ber);
	
	///////  serの計算と表示    /////////
	// ser=cs/(attempt*N);
	// putchar('\n');
	// printf(" ρ[dB]         SER \n");
	// printf(" %d          %.10f \n",rho,ser);
}

int main(void)
{	
	int i,j;

	////////////////      乱数の種を設定     ///////////////
	srand((unsigned int)time(NULL)); 
	
	/////////////// 時間と周波数の設定　　/////////////////////
	double t[N]; 						//時間t
	double f[N];						//周波数f
	
	for(i=0;i<N;i++){
		t[i]=T_sample*i;
		f[i]=delta_f*(i-N/2);
	}
	
	/////////   複素ベースバンド信号と搬送帯域信号の生成関数の試運転  /////
	double complex u_t[N];				//複素ベースバンド信号
	double complex s_t[N];				//搬送帯域信号
	
	// baseband_carrier(t,u_t,s_t);
	
	///////////////　　　　複素ガウス雑音の生成　　　　　　　////////////
	double complex awgn1[N];		
	
	boxmuller(awgn1,0,0.5);		//2σ^2=1より分散=電力=σ^2=0.5
	 
/* 	putchar('\n');
	printf("       ガウス雑音の表示 \n");
	printf("　時間　 　　    　実部　　　   　　 虚部 \n");
	for(i=0; i<N; i++){
	printf(" %f        " , t[i]);
	printf("%f          " , creal(awgn1[i]));
	printf("%f        \n" , cimag(awgn1[i]));
	}  */
	
	
	////////////   複素ガウス雑音の電力　　　　///////////////
	double P1t_awgn1=0.0;
	double P2t_awgn1=0.0;
	double P3t_awgn1=0.0;
	
/* 	 for(int k=0; k<loop; k++){	//loop=1と1000で行った
		boxmuller(awgn1,0,0.5);
		for(i=0; i<N; i++){
			// P1t_awgn1+=creal(awgn1[i])*creal(awgn1[i]);		//複素ガウス雑音のI軸の電力 ni^2
			// P2t_awgn1+=cimag(awgn1[i])*cimag(awgn1[i]);		//複素ガウス雑音のQ軸の電力 nq^2
			P3t_awgn1+=cabs(awgn1[i])*cabs(awgn1[i]);		//複素ガウス雑音の電力|n|^2　　n=ni+jnq 
		}
	}
	putchar('\n');
	// printf("%f   " , P1t_awgn1/(N*loop));
	// printf("%f   " , P2t_awgn1/(N*loop));
	printf("%f \n" , P3t_awgn1/(N*loop)); 
 */

	
	////////////　　複素ガウス雑音の確率密度関数   ///////////////
	int c_amp[loop],c_arg[loop];		//振幅，位相カウント用
	double amp[N],arg[N];				//振幅，位相格納用
	
	for(i=0; i<loop; i++){
		c_amp[i]=0;
		c_arg[i]=0;
	}
	
	
	int k;
	for(i=0; i<loop; i++){	
		for(j=0; j<N; j++){
			boxmuller(awgn1,0.0,0.5);		//(μ=0,σ^2=0.5)のガウス雑音
			amp[j]=cabs(awgn1[j]);
			arg[j]=atan2(creal(awgn1[j]),cimag(awgn1[j]))/M_PI;		
		}
		
		for(k=0; k<range; k++){				//振幅の分類 (0～0.05)(0.05～1.00)・・・()
			for(j=0; j<N; j++){
				if(k*w<=amp[j] && amp[j]<(k+1)*w) 
					c_amp[k]++;
			}
		}
		
		for(k=-range/2; k<range/2; k++){	//位相の分類
			for(j=0; j<N;j++){
				if(k*w<=arg[j] && arg[j]<(k+1)*w) 
					c_arg[k+range/2]++;
			}
		}
	}
	
	//////////　確率密度関数の表示  ///////////
	putchar('\n');
	printf("振幅の確率密度関数 \n");
	for(k=0; k<range; k++)
		printf(" %f    %f \n" , k*w,c_amp[k]/(N*loop*w));
	
	
	putchar('\n');
	printf("位相の確率密度関数 \n");	 
	for(k=-range/2; k<range/2; k++)	
		printf(" %f    %f \n" , k*w,c_arg[k+range/2]/(N*loop*w*M_PI));			//πで割る理由は横軸がradian/piだから全領域で積分して1になるようにπで割る．
	
	 
	 
	////////////　　　ガウス雑音のスペクトルの算出と表示　　　///////////////
	double complex fft_awgn1[N];			//ガウス雑音awgn1(μ=0,σ^2=0.5)のスペクトル
	double Pf_awgn1=0.0;					//周波数領域の電力
	
	/* for(int k=0; k<loop; k++){
	boxmuller(awgn1,0.0,0.5);
	FFT(awgn1,fft_awgn1);
	
 	// putchar('\n');
	// printf("　　　　ガウス雑音の振幅スペクトルを求める　　　\n");
	// printf("　周波数　　　    　振幅スペクトル\n");
	for(i=0; i<N; i++){
		//  周波数   //
		// printf("%f        " , f[i]);
		
		//　　ガウス雑音のスペクトル　　//
		// printf("%.5f   " , creal(fft_awgn1[i]));
		// printf("%.5f   " , cimag(fft_awgn1[i]));
		// printf("   %.17f      \n" , cabs(fft_awgn1[i]));
		Pf_awgn1+=cabs(fft_awgn1[i])*cabs(fft_awgn1[i]);
	} 
	}
	printf("%f \n",Pf_awgn1/(loop*N));   */
	
	/////////////　　受信信号とBER特性，SER特性　　　//////////////////
	double complex r_t[N];		//r(t)=s(t)+n(t)
	double complex awgn2[N];	//ガウス雑音awgn2(指導書後半の(3-2)式のσ^2から分散を計算)
	double sigma2;				//分散			
	int EbN0;					//ρ=Eb/N0
	
	putchar('\n');
	// printf(" ρ[dB]         BER \n");
	// for(EbN0=0; EbN0<=9; EbN0++){
	// 	sigma2=(double)BnTs/(2*s*c)*pow(10.0,(-0.1*EbN0));	
		
	// 	BER_QPSK(t,sigma2,EbN0);
	// }
	
	return 0;
}