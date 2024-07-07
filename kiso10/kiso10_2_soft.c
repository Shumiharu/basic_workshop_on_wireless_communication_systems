/*/////////////////////  NSC符号，RSC符号の軟判定復号(2乗ユークリッド距離を用いる)　　///////////////////////////////
・RSC符号の状態線図がNSC符号と全く違うことに注意
・
*/
#include <stdio.h>
#include <stdlib.h>		//乱数使用時のプログラム
#include <time.h>		//time関数(乱数の種を設定するため)使用時のプログラム
#include <math.h>		//数学的関数使用時のプログラム
#include <complex.h>	//複素数使用時のプログラム


#define a 1						//波数
#define s 1						//オーバーサンプル比
#define K 799			//シンボル，データ長
#define N (s*K)					//全サンプル数
#define sigmaN 0.5				//ガウス雑音の分散
#define Ts 0.001				//1シンボル時間
#define fc a/Ts 				//搬送波周波数fc
#define T_sample (Ts/s) 		//1サンプル当たりの時間
#define delta_f (1/(Ts*K))		//Δf=1/T Δf:スペクトルの幅  T(=Ts*K):時間長

//畳み込み符号用
#define delay_units 2			//遅延素子数
#define code_length 10			//符号長
#define r 2						//符号長を何倍にするか
#define code_rate 1.0/r 		//符号化率
#define choise_code 2				//1=NSC符号 2=RSC符号
//0=2乗ユークリッド距離，軟判定復号

//BER特性用↓↓↓↓
#define bit 2 				//1シンボル当たりのビット数 QPSKの場合は2
#define BnTs 1   				
#define first 50000
#define minEbN0 0
#define maxEbN0 12


////////////////////////関数プロトタイプ宣言/////////////////////////////////
void FFT(double complex Time[],double complex Freq[]);
void IFFT(double complex Freq[],double complex Time[]);
void QPSK(double complex base[]);
void inter(double complex s1[],double complex s2[]);		//先頭サンプルのみ抜き出す
void boxmuller(int num,double complex noise[],double mean,double variance);
void SP(double complex Serial[],double complex Parallel[]);
void PS(double complex Parallel[],double complex Serial[]);
void NSC(int nd,int input[],int output[]);
void RSC(int nd,int input[],int output[]);
void QPSK_bit(int nd,double complex input[],int output[]);
void bit_QPSK(int nd, int input[], double complex output[]);
double arg_min(double n1, double n2);		//☆☆☆☆ハミング距離ではint型だった
void RSC_state(int nd, int input[], int memory[]);



double arg_max(double n1, double n2){
	if(n1>=n2)
		return n1;		
	else
		return n2;
}

int main(void)
{
///////////////////////      乱数の種を設定     ///////////////
	srand((unsigned int)time(NULL)); 

	int i,j,m;
	double complex st[N];		//送信QPSK信号，データ数N
	int st_bit[N*bit];			//送信情報ビット，データ数N*bit
	int state[1];				//RSC符号の場合の最終状態の記録
	int st_tail[N*bit+2];		//送信情報ビットにテールビット(00)を挿入
	int codeword[N*bit*r+2*r];		//符号語，要素がr倍される
	double complex codeword_QPSK[N*r+r];	//符号語の複素QPSK信号[N*bit*r/bit+2*r/bit=N*r+r]
	double complex awgn[N*r+r];			//AWGN
	double complex rt_QPSK[N*r+r];		//rt_QPSK=AWGN+codeword_QPSK
	int rt_bit[N*bit*r+2*r];				//rt_bit=Bit(rt_QPSK)
	double pass_metA[2][N*r+r];		//時刻tにおける状態Aでのパスメトリック
	double pass_metB[2][N*r+r];		//時刻tにおける状態Bでのパスメトリック
	double pass_metC[2][N*r+r];		//時刻tにおける状態Cでのパスメトリック
	double pass_metD[2][N*r+r];		//時刻tにおける状態Dでのパスメトリック
	int key;						//0=A(00)，1=B(10)，2=C(01)，3=D(11)
	int estimate_tail[N*bit+2];
	int estimate_bit[N*bit];			//推定送信情報ビット系列
	double complex estimate_QPSK[N];	//推定QPSK複素信号
	
	double amp=1/(sqrt(2*s));
	double complex QPSK_signal[4]={amp+I*amp,-amp+I*amp,-amp-I*amp,amp-I*amp};
	
///////////////////////////////////BER特性算出時に用いる変数　　/////////////////////////////////////////////////
	double sigma;				//ガウス雑音の分散			
	int EbN0;					//Eb/N0
	double ci,cq;				//BERカウント用
	double ber;  				//誤り率
	int loop,LOOP;
	
	
	if(choise_code==1)	printf("NSC符号\n");
	if(choise_code==2)	printf("RSC符号\n");
	printf("符号長=%d, 符号化率=%f, 遅延素子数=%d, s=%d\n",r*(N*bit+delay_units),code_rate,delay_units,s);
	printf(" ρ[dB]         BER \n");
	for(EbN0=minEbN0; EbN0<=maxEbN0; EbN0++){
		ci=0.0; cq=0.0; 
		LOOP=first;			//Eb/N0の値ごとにfirst回でBERを軽く算出する
		//☆☆☆☆☆☆☆☆☆雑音の大きさ
		sigma=(double)BnTs/(2*s*bit)*pow(10.0,(-0.1*EbN0))*((2*N*bit+2*delay_units)/(N*bit));		//ガウス雑音の分散算出
		for(loop=0; loop<LOOP; loop++){	
			
			QPSK(st);
			QPSK_bit(N,st,st_bit);
			
			//ビットをコピー
			for(i=0; i<N*bit; i++)
				st_tail[i]=st_bit[i];
			
			//テールビット挿入(NSC符号)
			if(choise_code==1){		//NSC符号の場合は最終状態がどこであってもテールビットとして(00)を挿入すればAに終端する
				st_tail[N*bit]=0;	
				st_tail[N*bit+1]=0;
			}
			
			//テールビット挿入(RSC符号)
			if(choise_code==2){
				RSC_state(N*bit,st_bit,state);		//RSC符号の場合は最終状態によってテールビットが異なるから，まず情報ビットだけで作られる符号語から最終状態を記録する
		
				if(state[0]==0){	//状態Aにいるとき以下のテールビットを挿入し，状態Aで終端させる
					st_tail[N*bit]=0;	
					st_tail[N*bit+1]=0;
				}
				
				if(state[0]==1){	//状態Bにいるとき以下のテールビットを挿入し，状態Aで終端させる
					st_tail[N*bit]=1;	
					st_tail[N*bit+1]=1;
				}
				
				if(state[0]==2){	//状態Cにいるとき以下のテールビットを挿入し，状態Aで終端させる
					st_tail[N*bit]=1;	
					st_tail[N*bit+1]=0;
				}
				
				if(state[0]==3){	//状態Dにいるとき以下のテールビットを挿入し，状態Aで終端させる
					st_tail[N*bit]=0;	
					st_tail[N*bit+1]=1;
				}
			}
			
			//符号の選択
			if(choise_code==1)	NSC(N*bit+2,st_tail,codeword);
			if(choise_code==2)	RSC(N*bit+2,st_tail,codeword);
			
			bit_QPSK(N*bit*r+2*r,codeword,codeword_QPSK);
			
			///////　ガウス雑音付加　　////////////////////////////////////////////////////////////////////////
			boxmuller(N*r+r,awgn,0.0,sigma);		//　AWGNの生成　
			
			for(i=0; i<N*r+r; i++){		//☆☆☆雑音のくわえ方＿テールビットにも雑音くわえる
				rt_QPSK[i]=codeword_QPSK[i]+awgn[i];		
			}
			
			//2乗ユークリッド距離を用いる場合は入力するのは複素数
			QPSK_bit(N*r+r,rt_QPSK,rt_bit);
			
			//パスメトリック計算
			for(j=0; j<2; j++){
				for(i=0; i<N*r+r; i++){
					pass_metA[j][i]=0.0;
					pass_metB[j][i]=0.0;
					pass_metC[j][i]=0.0;
					pass_metD[j][i]=0.0;
				}
			}
			
			//NSC符号の場合
			if(choise_code==1){
				//軟判定ビタビ復号_パスメトリック計算
				for(i=0; i<N*r+r; i++){
					if(i==0){
						for(j=0; j<2; j++){
							pass_metA[j][i]=cabs(rt_QPSK[i]-QPSK_signal[0]);	//A→A (00)=(A+jA)とのユークリッド距離
							pass_metB[j][i]=cabs(rt_QPSK[i]-QPSK_signal[2]);	//A→B (11)=(-A-jA)とのユークリッド距離
						}
					}
					if(i==1){
						for(j=0; j<2; j++){
							pass_metA[j][i]=pass_metA[0][i-1]+cabs(rt_QPSK[i]-QPSK_signal[0]);	//A→A passA+(00)=(A+jA)とのユークリッド距離
							pass_metB[j][i]=pass_metA[0][i-1]+cabs(rt_QPSK[i]-QPSK_signal[2]);	//A→B passA+(11)=(-A-jA)とのユークリッド距離
							pass_metC[j][i]=pass_metB[0][i-1]+cabs(rt_QPSK[i]-QPSK_signal[3]);	//B→C passB+(01)=(A-jA)とのユークリッド距離
							pass_metD[j][i]=pass_metB[0][i-1]+cabs(rt_QPSK[i]-QPSK_signal[1]);	//B→D passB+(10)=(-A+jA)とのユークリッド距離	
						}
					}
					if(i>=2){					
						pass_metA[0][i]=arg_min(pass_metA[0][i-1],pass_metA[1][i-1])+cabs(rt_QPSK[i]-QPSK_signal[0]);	//A→A passA+(00)=(A+jA)とのユークリッド距離
						pass_metA[1][i]=arg_min(pass_metC[0][i-1],pass_metC[1][i-1])+cabs(rt_QPSK[i]-QPSK_signal[2]);	//C→A passC+(11)=(-A-jA)とのユークリッド距離
						pass_metB[0][i]=arg_min(pass_metA[0][i-1],pass_metA[1][i-1])+cabs(rt_QPSK[i]-QPSK_signal[2]);	//A→B passA+(11)=(-A-jA)とのユークリッド距離
						pass_metB[1][i]=arg_min(pass_metC[0][i-1],pass_metC[1][i-1])+cabs(rt_QPSK[i]-QPSK_signal[0]);	//C→B passC+(00)=(A+jA)とのユークリッド距離
						pass_metC[0][i]=arg_min(pass_metB[0][i-1],pass_metB[1][i-1])+cabs(rt_QPSK[i]-QPSK_signal[3]);	//B→C passB+(01)=(A-jA)とのユークリッド距離
						pass_metC[1][i]=arg_min(pass_metD[0][i-1],pass_metD[1][i-1])+cabs(rt_QPSK[i]-QPSK_signal[1]);	//D→C passD(10)=(-A+jA)とのユークリッド距離
						pass_metD[0][i]=arg_min(pass_metB[0][i-1],pass_metB[1][i-1])+cabs(rt_QPSK[i]-QPSK_signal[1]);	//B→D passB+(10)=(-A+jA)とのユークリッド距離
						pass_metD[1][i]=arg_min(pass_metD[0][i-1],pass_metD[1][i-1])+cabs(rt_QPSK[i]-QPSK_signal[3]);	//D→D passD+(01)=(A-jA)とのユークリッド距離	
					}
					
				}
				
				//軟判定ビタビ復号_ビタビアルゴリズム
				key=0;		//最初はAから判定	
				//2.2乗ユークリッド距離が最大となる状態からさかのぼっていき，送信情報を推定
				for(i=(N*bit+2)-1; i>=0; i--){
					if(key==0){		//状態A(00)での判定
						if(pass_metA[0][i]<=pass_metA[1][i])	key=0;	//次回状態Aで判定(=をこっちにつけることによって時刻t=0,1もうまくいく気がする)
						if(pass_metA[0][i]>pass_metA[1][i])		key=2;	//次回状態Cで判定
						estimate_tail[i]=0;		//状態Aでの判定は入力が必ず0
					}else if(key==1){		//状態B(10)での判定
						if(pass_metB[0][i]<=pass_metB[1][i])	key=0;	//次回状態Aで判定
						if(pass_metB[0][i]>pass_metB[1][i])		key=2;	//次回状態Cで判定
						estimate_tail[i]=1;		//状態Bでの判定は入力が必ず1
					}else if(key==2){		//状態C(01)での判定
						if(pass_metC[0][i]<=pass_metC[1][i])	key=1;	//次回状態Bで判定
						if(pass_metC[0][i]>pass_metC[1][i])		key=3;	//次回状態Dで判定
						estimate_tail[i]=0;		//状態Cでの判定は入力が必ず0
					}else if(key==3){		//状態D(11)での判定
						if(pass_metD[0][i]<=pass_metD[1][i])	key=1;	//次回状態Bで判定
						if(pass_metD[0][i]>pass_metD[1][i])		key=3;	//次回状態Dで判定
						estimate_tail[i]=1;		//状態Bでの判定は入力が必ず1
					}	
					// printf(" key=%d  ",key);	
				}
			}//NSC符号_ビタビ復号 end
					
			
			
			//RSC符号の場合
			if(choise_code==2){
				//軟判定ビタビ復号_パスメトリック計算
				for(i=0; i<N*r+r; i++){
					if(i==0){
						for(j=0; j<2; j++){
							pass_metA[j][i]=cabs(rt_QPSK[i]-QPSK_signal[0]);	//A→A (00)=(A+jA)とのユークリッド距離
							pass_metB[j][i]=cabs(rt_QPSK[i]-QPSK_signal[2]);	//A→B (11)=(-A-jA)とのユークリッド距離
						}
					}
					if(i==1){
						for(j=0; j<2; j++){
							pass_metA[j][i]=pass_metA[0][i-1]+cabs(rt_QPSK[i]-QPSK_signal[0]);	//A→A passA+(00)=(A+jA)とのユークリッド距離
							pass_metB[j][i]=pass_metA[0][i-1]+cabs(rt_QPSK[i]-QPSK_signal[2]);	//A→B passA+(11)=(-A-jA)とのユークリッド距離
							pass_metC[j][i]=pass_metB[0][i-1]+cabs(rt_QPSK[i]-QPSK_signal[1]);	//B→C passB+(10)=(-A+jA)とのユークリッド距離
							pass_metD[j][i]=pass_metB[0][i-1]+cabs(rt_QPSK[i]-QPSK_signal[3]);	//B→D passB+(01)=(A-jA)とのユークリッド距離	
						}
					}
					if(i>=2){					
						pass_metA[0][i]=arg_min(pass_metA[0][i-1],pass_metA[1][i-1])+cabs(rt_QPSK[i]-QPSK_signal[0]);	//A→A passA+(00)=(A+jA)とのユークリッド距離
						pass_metA[1][i]=arg_min(pass_metC[0][i-1],pass_metC[1][i-1])+cabs(rt_QPSK[i]-QPSK_signal[2]);	//C→A passC+(11)=(-A-jA)とのユークリッド距離
						pass_metB[0][i]=arg_min(pass_metA[0][i-1],pass_metA[1][i-1])+cabs(rt_QPSK[i]-QPSK_signal[2]);	//A→B passA+(11)=(-A-jA)とのユークリッド距離
						pass_metB[1][i]=arg_min(pass_metC[0][i-1],pass_metC[1][i-1])+cabs(rt_QPSK[i]-QPSK_signal[0]);	//C→B passC+(00)=(A+jA)とのユークリッド距離
						pass_metC[0][i]=arg_min(pass_metB[0][i-1],pass_metB[1][i-1])+cabs(rt_QPSK[i]-QPSK_signal[1]);	//B→C passB+(10)=(-A+jA)とのユークリッド距離
						pass_metC[1][i]=arg_min(pass_metD[0][i-1],pass_metD[1][i-1])+cabs(rt_QPSK[i]-QPSK_signal[3]);	//D→C passD+(01)=(A-jA)とのユークリッド距離
						pass_metD[0][i]=arg_min(pass_metB[0][i-1],pass_metB[1][i-1])+cabs(rt_QPSK[i]-QPSK_signal[3]);	//B→D passB+(01)=(A-jA)とのユークリッド距離
						pass_metD[1][i]=arg_min(pass_metD[0][i-1],pass_metD[1][i-1])+cabs(rt_QPSK[i]-QPSK_signal[1]);	//D→D passD+(10)=(-A+jA)とのユークリッド距離	
					}
				}
				
				//軟判定ビタビ復号_ビタビアルゴリズム	
				key=0;		//最初は状態Aから判定
				for(i=(N*bit+2)-1; i>=0; i--){
					if(key==0){		//状態A(00)での判定
						if(pass_metA[0][i]<=pass_metA[1][i]){	
							key=0;	//次回状態Aで判定(=をこっちにつけることによって時刻t=0,1もうまくいく気がする)
							estimate_tail[i]=0;		//A→Aは入力0
						}
						if(pass_metA[0][i]>pass_metA[1][i]){
							key=2;	//次回状態Cで判定
							estimate_tail[i]=1;		//C→Aは入力1
						}
					}else if(key==1){		//状態B(10)での判定
						if(pass_metB[0][i]<=pass_metB[1][i]){
							key=0;	//次回状態Aで判定
							estimate_tail[i]=1;		//A→Bは入力1
						}
						if(pass_metB[0][i]>pass_metB[1][i]){
							key=2;	//次回状態Cで判定
							estimate_tail[i]=0;		//C→Bは入力0
						}
					}else if(key==2){		//状態C(01)での判定
						if(pass_metC[0][i]<=pass_metC[1][i]){
							key=1;	//次回状態Bで判定
							estimate_tail[i]=1;		//B→Cは入力1
						}
						if(pass_metC[0][i]>pass_metC[1][i]){
							key=3;	//次回状態Dで判定
							estimate_tail[i]=0;		//D→Cは入力0
						}	
					}else if(key==3){		//状態D(11)での判定
						if(pass_metD[0][i]<=pass_metD[1][i]){
							key=1;	//次回状態Bで判定
							estimate_tail[i]=0;		//B→Dは入力0
						}
						if(pass_metD[0][i]>pass_metD[1][i]){
							key=3;	//次回状態Dで判定
							estimate_tail[i]=1;		//D→Dは入力1
						}
					}	
				}
			}//RSC符号_ビタビ復号 end

			 
			 
			// tail bit removal
			for(i=0; i<N*bit; i++)
				estimate_bit[i]=estimate_tail[i];


			//ビット系列から複素信号
			bit_QPSK(N*bit,estimate_bit,estimate_QPSK);
			
			/// 誤り判定(先頭点のみ)  ////////////////////////////////////////////////////////////////////////
			for(j=0; j<K; j++){
				if(creal(st[j*s])>=0 && creal(estimate_QPSK[j*s])<0) ci++;
				if(creal(st[j*s])<0 && creal(estimate_QPSK[j*s])>=0) ci++;
				if(cimag(st[j*s])>=0 && cimag(estimate_QPSK[j*s])<0) cq++;	
				if(cimag(st[j*s])<0 && cimag(estimate_QPSK[j*s])>=0) cq++;
			}	
				//BERの値によってループ調整
			if(loop==LOOP-1){
				ber=(ci+cq)/(bit*(double)(1+loop)*K);		//一旦BERを算出
				
				if(ber>=1E-1 && ber<1)
					LOOP=pow(10,4);			//loop=first=5000→このifに入る→LOOP=10が代入される→最初のforの部分でloop++される→loop<LOOPの判定でloopの方が大きいのでこのforを抜け出せる→BERを計算して結果表示
				if(ber>=1E-2 && ber<1E-1)
					LOOP=pow(10,4);
				if(ber>=1E-3 && ber<1E-2)
					LOOP=pow(10,4);
				if(ber>=1E-4 && ber<1E-3)
					LOOP=5*pow(10,5);
				if(ber>=1E-5 && ber<1E-4)
					LOOP=5*pow(10,6);
				if(ber>=1E-6 && ber<1E-5)
					LOOP=pow(10,7);
				if(ber>=1E-7 && ber<1E-6)
					LOOP=5*pow(10,7);
				if(ber>=1E-8 && ber<1E-7)
					LOOP=5*pow(10,7);
				if(ber>=0 && ber<1E-8)
					LOOP=pow(10,8);
				// printf("BER=%e, loop=%d, LOOP=%d, Eb/N0=%d dB\n",ber,loop,LOOP,EbN0);
			}

			ber=(ci+cq)/(bit*(double)(1+loop)*K);
			printf("loop=%d         Eb/N0=%d         %e\r",loop,EbN0,ber);	
			
		
			
		}//loop end
			
		ber=(ci+cq)/(bit*(double)loop*K);
		printf("loop=%d         Eb/N0=%d         %e\n",loop,EbN0,ber);	
			
			
	}//Eb/N0 end
	
	
	return 0;
}

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

////////////  高速フーリエ逆変換IFFT(Inverse FFT)　　/////////////
void IFFT(double complex Freq[],double complex Time[]){
	int i,j,k,m;
	double arg;
	double complex X[N],Y[N];
	
	for(i=0;i<N/2;i++){
		X[i]=Freq[i+N/2];
		X[i+N/2]=Freq[i];
		Time[i]=0;
		Time[i+N/2];
	}
	
	
	for(m=2;m<=N;m*=2){
		j=m/2;
		for(i=0;i<N/m;i++){
			for(k=0;k<j;k++){
				arg=2*M_PI*(k)/m;
				Y[i*m+k]=X[i*j+k]+cexp(arg*I)*X[i*j+(N/2)+k];//バタフライ上側
				Y[i*m+j+k]=X[i*j+k]-cexp(arg*I)*X[i*j+(N/2)+k];//バタフライ下側
			}
		}
		for(i=0;i<N;i++){
			X[i]=Y[i];
		}
		for(i=0;i<N;i++){
			Y[i]=0.0+0.0*I;
		}
	}
	for(i=0;i<N;i++){
		Time[i]=X[i]/sqrt(N);
	}
}

/////////////   複素ベースバンド信号(と搬送帯域信号)の生成　　　//////////////////////
void QPSK(double complex base[]){
		int i,j;
		double tmp_I;	//同相成分Iを一時的に格納
		double tmp_Q;	//直交成分Qを一時的に格納
		
		for(i=0; i<K; i++){
			tmp_I=1/(sqrt(2*s))*(-2*(rand()%2)+1);	//電圧の大きさが1/√2s , 乱数が0なら１で1ならー1	
			tmp_Q=1/(sqrt(2*s))*(-2*(rand()%2)+1);	//電圧の大きさが1/√2s , 乱数が0なら１で1ならー1
			
			for(j=0; j<s; j++)
				base[s*i+j]=tmp_I+tmp_Q*I;				//複素ベースバンド信号u(t)=I(t)+jQ(t)
		}			
}

////////////////　先頭サンプルのみ抜き出し関数　/////////////////////////
void inter(double complex s1[],double complex s2[]){
	int i,j;
	
	for(i=0;i<N;i+=s){		//シンボル数だけloop
			s2[i]=s1[i];		
			for(j=1;j<s;j++)	//1サンプル以降を0にする
				s2[i+j]=0.0+0.0*I;
		}
}

//　　　　Box-Muller法によるnum個からなる配列=ガウス雑音生成　　　　　　//
void boxmuller(int num,double complex noise[], double mean , double variance){
	int i,j;
	double u1,u2;
	double nreal[num],nimag[num];
	
	for(i=0; i<num; i++){
	////////////////　　　　[0,1]の一様乱数生成　　　　　/////////////

		u1=((double)rand()+1.0)/((double)RAND_MAX+1.0);		//rand()=0～32767
		u2=((double)rand()+1.0)/((double)RAND_MAX+1.0);
		
		///   標準偏差σを1/√2にすれば各軸の分散σ^2=1/2になる．
		nreal[i]=sqrt(variance)*sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2)+mean;		//分散を1/√2にして合成後の電力を1Wにする．
		nimag[i]=sqrt(variance)*sqrt(-2.0*log(u1))*sin(2.0*M_PI*u2)+mean;
		noise[i]=nreal[i]+nimag[i]*I;
	}
}
/////////////////////////Serial to Pallel/////////////////////////
void SP(double complex Serial[],double complex Parallel[]){
	int i,j;
	
	j=0;
	for(i=0; i<N; i++){
		if(i>=(N-K)/2 && i<(N+K)/2){
			Parallel[i]=Serial[j*s];
			j++;
		}
		else{
			Parallel[i]=0.0+0.0*I;
		}
	}	
}

////////////////////////Parallel to Serial///////////////////////
void PS(double complex Parallel[],double complex Serial[]){
	int i,j;
	
	j=0;
	for(i=0; i<N; i++){
		if(i>=(N-K)/2 && i<(N+K)/2){
			Serial[j*s]=Parallel[i];
			j++;
		}
	}	
}

///////////////////　畳み込み符号NSC符号生成 /////////////////////////
//                入力[nd]	出力[符号化率^-1*nd]
void NSC(int nd,int input[],int output[]){
	int i,j,k,n1,n2;
	int D[delay_units][nd];			//D[0][]=D0，D[1][]=D1，D[2][]=D2，...
	int v1[nd],v2[nd];				
	
	
	//遅延素子の初期化
	for(i=0;i<delay_units;i++)
		D[i][0]=0;		//初期値として0を入れる
	
	/* //遅延素子に値を入れる
	for(i=0;i<delay_units;i++){
		for(j=1;j<nd;j++){		//遅延素子の先頭には0が入ってるからj=1からスタート	
			if(i==0) 
				D[i][j]=input[j-1];		//入力系列の1時刻前のデータを記憶させる
			else
				D[i][j]=D[i-1][j-1];	//1つ前の遅延素子の状態を最初から入れる
		}
	} */
	
	//遅延素子に値を入れる
	for(j=1;j<nd;j++){		//遅延素子の先頭には0が入ってるからj=1からスタート	
		for(i=0;i<delay_units;i++){
			if(i==0) 
				D[i][j]=input[j-1];		//入力系列の1時刻前のデータを記憶させる
			else
				D[i][j]=D[i-1][j-1];	//1つ前の遅延素子の状態を最初から入れる
		}
	}
	//NSC符号器の2つの出力を計算
	for(j=0;j<nd;j++){
		v1[j]=(input[j]+D[1][j])%2;				//V1=U+D2のmod2(=XOR)
		v2[j]=(input[j]+D[0][j]+D[1][j])%2;		//V2=U+D1+D2のmod2(=XOR)
	}
	
	//符号語outputを作成
	n1=0; n2=0;
	for(k=0;k<nd*r;k++){
			if(k%2==0){	
				output[k]=v1[n1];	//kが偶数番目はv1を
				n1++;
			}
			if(k%2==1){
				output[k]=v2[n2];	//kガキ数番目はv2を
				n2++;
			}
	}
	
	//遅延素子の状態確認
	// for(j=0;j<delay_units;j++){
		// putchar('\n');
		// for(i=0;i<nd;i++)
			// printf("%d",D[j][i]);
	// }
	// putchar('\n');
}

///////////////////　畳み込み符号RSC符号生成 /////////////////////////
//				   入力[nd]	出力[符号化率/1*nd]
void RSC(int nd,int input[],int output[]){
	int i,j,k,n1,n2;
	int D[delay_units][nd];			//D[0][]=D0，D[1][]=D1，D[2][]=D2，...
	int v1[nd],v2[nd];				
	
	
	//遅延素子の初期化
	for(i=0;i<delay_units;i++)
		D[i][0]=0;		//初期値として0を入れる
	
	/* //遅延素子に値を入れる  
	for(i=0;i<delay_units;i++){
		for(j=1;j<nd;j++){		//遅延素子の先頭には0が入ってるからj=1からスタート	
			if(i==0){ 
				D[i][j]=(input[j-1]+D[0][j-1]+D[1][j-1])%2;		//入力系列の1時刻前のデータを記憶させる
				                              ----------
											  この部分が最初i=0の遅延素子の状態を計算するからまだ値が入ってない状態になっていた
											  →forの並びを変える
			}
			else{
				D[i][j]=D[i-1][j-1];	//D2[t]=D1[t-1]，1つ前の遅延素子の状態を最初から入れる
	
			}
		}
	} */
	
	//遅延素子に値を入れる
	
	for(j=1;j<nd;j++){		//遅延素子の先頭には0が入ってるからj=1からスタート	
		for(i=0;i<delay_units;i++){		
			if(i==0){ 
				D[i][j]=(input[j-1]+D[0][j-1]+D[1][j-1])%2;		//入力系列の1時刻前のデータを記憶させる
			}
			else{
				D[i][j]=D[i-1][j-1];	//D2[t]=D1[t-1]，1つ前の遅延素子の状態を最初から入れる
			}
		}
	}
	
	//RSC符号器の2つの出力を計算
	for(j=0;j<nd;j++){
		v1[j]=input[j];							///V1=U（情報ビットそのまま）
		v2[j]=(input[j]+D[0][j]+D[1][j]+D[1][j])%2;		//V2=U+D1+D2+D2のmod2(=XOR)
	}
	
	//符号語outputを作成
	n1=0; n2=0;
	for(k=0;k<nd*r;k++){
		if(k%2==0){	
			output[k]=v1[n1];	//kが偶数番目はv1を
			n1++;
		}
		if(k%2==1){
			output[k]=v2[n2];	//kガキ数番目はv2を
			n2++;
		}
	}
	
	//遅延素子の状態の確認
	// for(j=0;j<delay_units;j++){
		// putchar('\n');
		// for(i=0;i<nd;i++)
			// printf("%d",D[j][i]);
	// }
}

/////////////////// QPSK複素信号をビット系列を変換する関数 //////////////////
//             入力のデータ数     QPSK複素信号(nd)   ビット系列(nd*bit)
void QPSK_bit(int nd,double complex input[],int output[]){
	int i,j;
	
	for(i=0;i<nd;i++){
		if(creal(input[i])>0 && cimag(input[i])>0){		//第1象限(00)
			output[i*bit]=0;
			output[i*bit+1]=0;
		}
		if(creal(input[i])<0 && cimag(input[i])>0){		//第2象限(10)
			output[i*bit]=1;
			output[i*bit+1]=0;
		}
		if(creal(input[i])<0 && cimag(input[i])<0){		//第3象限(11)
			output[i*bit]=1;
			output[i*bit+1]=1;
		}
		if(creal(input[i])>0 && cimag(input[i])<0){		//第4象限(01)
			output[i*bit]=0;
			output[i*bit+1]=1;
		}	
	}
}

/////////////////////　ビット列をQPSK信号に変換する関数  //////////////////////
//			    ビット長   ビット系列(nd)        QPSK信号(nd/bit) 
void bit_QPSK(int nd, int input[], double complex output[]){
	int i,j;
	
	for(i=0;i<nd/bit; i++){
		if(input[i*bit]==0 && input[i*bit+1]==0)
			output[i]=1/(sqrt(2*s))+I*1/(sqrt(2*s));	//第1象限(00)=A+jA
		
		if(input[i*bit]==1 && input[i*bit+1]==0)
			output[i]=-1/(sqrt(2*s))+I*1/(sqrt(2*s));	//第2象限(10)=-A+jA
		
		if(input[i*bit]==1 && input[i*bit+1]==1)
			output[i]=-1/(sqrt(2*s))-I*1/(sqrt(2*s));	//第1象限(11)=-A-jA
		
		if(input[i*bit]==0 && input[i*bit+1]==1)
			output[i]=1/(sqrt(2*s))-I*1/(sqrt(2*s));	//第1象限(01)=A-jA
	}
}


double arg_min(double n1, double n2){
	if(n1<=n2)
		return n1;		
	else
		return n2;
}
	
////////////////RSC符号の最終状態を記憶させる関数//////////////////////////////
//               　　ビット長            情報ビット      状態を記録
void RSC_state(int nd, int input[], int memory[]){
	// int hugougo[data_length*r];		//符号語,ビット数がr倍される
	
	
	int i,j,k,n1,n2;
	int D[delay_units][nd];			//D[0][]=D0，D[1][]=D1，D[2][]=D2，...
	int v1[nd],v2[nd];				
	
	
	//遅延素子の初期化
	for(i=0;i<delay_units;i++)
		D[i][0]=0;		//初期値として0を入れる
	
	
	//遅延素子に値を入れる
	for(j=1; j<nd ;j++){		//遅延素子の先頭には0が入ってるからj=1からスタート	
		for(i=0; i<delay_units; i++){		
			if(i==0) 
				D[i][j]=(input[j-1]+D[0][j-1]+D[1][j-1])%2;		//入力系列の1時刻前のデータを記憶させる
			else
				D[i][j]=D[i-1][j-1];	//D2[t]=D1[t-1]，1つ前の遅延素子の状態を最初から入れる		
		}
	}
	
	//最終時刻の遅延素子の状態から最終的な状態をmemoryに記憶
	if(D[0][nd-1]==0 && D[1][nd-1]==0){		//終端手前の状態がAのとき
		if(input[nd-1]==0)	memory[0]=0;	//最終状態はA
		if(input[nd-1]==1)	memory[0]=1;	//最終状態はB
	}
	
	if(D[0][nd-1]==1 && D[1][nd-1]==0){		//終端手前の状態がBのとき
		if(input[nd-1]==0)	memory[0]=3;	//最終状態はD
		if(input[nd-1]==1)	memory[0]=2;	//最終状態はC
	}
	
	if(D[0][nd-1]==0 && D[1][nd-1]==1){		//終端手前の状態がCのとき
		if(input[nd-1]==0)	memory[0]=1;	//最終状態はB
		if(input[nd-1]==1)	memory[0]=0;	//最終状態はA
	}
	
	if(D[0][nd-1]==1 && D[1][nd-1]==1){		//終端手前の状態がDのとき
		if(input[nd-1]==0)	memory[0]=2;	//最終状態はC
		if(input[nd-1]==1)	memory[0]=3;	//最終状態はD
	}
	//遅延素子の状態の確認
	// for(j=0;j<delay_units;j++){
		// putchar('\n');
		// for(i=0;i<nd;i++)
			// printf("%d",D[j][i]);
	// }
	
	//最終状態を確認
	// putchar('\n');
	// printf("情報ビットの最後のビット:%d\t最終状態：%d\n",input[nd-1],memory[0]);
}


	










