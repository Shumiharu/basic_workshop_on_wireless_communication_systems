////////////    関数化後のOFDM-QPSK伝送のBER特性     ///////////

#include <stdio.h>
#include <stdlib.h>		//乱数使用時のプログラム
#include <time.h>		//time関数(乱数の種を設定するため)使用時のプログラム
#include <math.h>		//数学的関数使用時のプログラム
#include <complex.h>	//複素数使用時のプログラム


#define a 1			 			//波数
#define s 1						//オーバーサンプル比
#define K 256					//サブキャリア数
#define N (s*K)					//全サンプル数
#define sigmaN 0.5				//ガウス雑音の分散
#define Ts 0.001				//1サンプル時間(OFDMでは)
#define fc a/Ts 				//搬送波周波数fc
#define T_sample (Ts/s) 		//1サンプル当たりの時間
#define delta_f (1/(Ts*K))		//Δf=1/T Δf:スペクトルの幅  T(=Ts*K):時間長

//OFDM特有のパラメータ↓↓↓↓
#define L 16						//L波マルチパス
#define dB 1					//何dB減衰
#define GI 32					//ガードインターバルのサンプル数


//BER特性用↓↓↓↓
#define bit 2.0 				//1シンボル当たりのビット数 QPSKの場合は2
#define BnTs 1   				
#define frame pow(10,3)			//フレーム数

////////////////////////関数プロトタイプ宣言/////////////////////////////////
void FFT(double complex Time[],double complex Freq[]);
void IFFT(double complex Freq[],double complex Time[]);
void QPSK(double complex base[]);
void QAM16(double complex QAM[]);
void inter(double complex s1[],double complex s2[]);		//先頭サンプルのみ抜き出す
void boxmuller(int num,double complex noise[],double mean,double variance);
void SP(double complex Serial[],double complex Parallel[]);
void PS(double complex Parallel[],double complex Serial[]);
void GI_insertion(double complex data[],double complex data_GI[]);
void GI_removal(double complex data_GI[],double complex data[]);
double average2(double complex data[]);

int main(void)
{	
	int i,j,m; 
	
///////////////////////      乱数の種を設定     ///////////////
	// srand((unsigned int)time(NULL)); 

////////////////////// 時間と周波数の設定　　//////////////////////
	double t[N]; 						//時間t    0～K*Tsまでの時間
	double f[N];						//周波数f  -Ts/s～Ts/sまでの周波数
	
	for(i=0;i<N;i++){
		t[i]=T_sample*i;
		f[i]=delta_f*(i-N/2);
	}
	
	
///////////////////// 準静的dB減衰L波マルチパスフェージングモデルの作成  ////////////////////////
	double complex ht_multifd[N];		//マルチパスフェージングのインパルス応答
	double complex Hf_multifd[N];		//マルチパスフェージングの伝達関数
	double P_multifd[L];				//各パスの電力
	double complex awgn1[1];			
	double Psum1=0.0;					//Hf_multifdの総電力
	
	double complex sF1[N];		//QPSK信号
	double complex sF2[N];		//先頭サンプルのみQPSK信号
	double complex sF[N];		//S/P変換後の信号
	double complex st[N];		//IFFT後の時間信号（雑音状）	
	
	double complex st_GI[(K+GI)*s];	//GI挿入後の時間信号
	double complex sF_GI[(K+GI)*s];	//GI挿入後の周波数領域信号
	double complex st3[(K+GI)*s];	//マルチパス通信路通過後の時間信号
	double complex Delay[(K+GI)*s];			//
	double complex tmp[N];				//
	double complex awgn[(K+GI)*s];		//ガウス雑音
	double complex rt1[(K+GI)*s];		//ガウス雑音重畳後の時間信号
	double complex rt2[N];				//GI除去後の時間信号
	double complex rF2[N];				//rt2のフーリエ変換
	double complex rhatF[N];			//FDE後の周波数信号
	double complex rF[N];				//rhatをP/S変換
	
	
///////////////////////////////////BER特性算出時に用いる変数　　/////////////////////////////////////////////////
	double sigma;				//ガウス雑音の分散			
	int EbN0;					//Eb/N0
	double ci,cq;				//BERカウント用
	double ber;  				//誤り率
	
	printf("%dパス準静的レイリーフェージング，%ddB減衰モデル　s=%d\n",L,dB,s);
	printf(" ρ[dB]         BER \n");
	for(EbN0=0; EbN0<=30; EbN0+=5){
		ci=0.0; cq=0.0; 
		sigma=(double)BnTs/(2*s*bit)*pow(10.0,(-0.1*EbN0))*((K+GI)/(double)K);		//GIの挿入によって帯域幅がK+GI/K倍拡大→雑音電力K+GI/K倍される
		//  遅延波格納用tmp　初期化 /////////////////////////////////////////////////////////////　　
		for(i=0;i<N;i++)
			tmp[i]=0+0*I;
		
		/// 　フレーム数分OFDMを伝送する　///////////////////////////////////////////////////
		for(m=0; m<frame; m++){	
			//// フレームごとに送信信号stを作成 //////////////////////////////////////////////////////////////////
			QPSK(sF1);			//QPSK信号生成
			inter(sF1,sF2);		//K本のサブキャリア=先頭サンプルのみ抜き出す
			SP(sF2,sF);			//
			IFFT(sF,st);
			//マルチパスフェージング作成
			
			P_multifd[0]=(1-pow(10,-(dB/10.0)))/(1-pow(10,-(dB*L/10.0)));
			for(i=0;i<K;i++){
				for(j=0;j<s;j++){				
					if(i<L && j==0){		//L波分作成　先頭にのみインパルスの値を入れる
						P_multifd[i]=P_multifd[0]*pow(10,-(dB*i/10.0));		//電力を減衰させてゆく P_multifd[0]=P_multifd[0]*1となる．
						boxmuller(1,awgn1,0.0,P_multifd[i]/2.0);			//2σ^2=P_multifd[i]となるようにガウス雑音を作成	
						ht_multifd[j+s*i]=awgn1[0];
					}else
						ht_multifd[j+s*i]=0+0*I;
				}
			}
	
			FFT(ht_multifd,Hf_multifd);		//マルチパスフェージングの伝達関数H(f)
			
			//H(f)の総電力計算→正規化の準備
			Psum1=0.0;		//☆☆☆
			for(i=0;i<N;i++)
				Psum1+=cabs(Hf_multifd[i])*cabs(Hf_multifd[i]);
			
			//平均電力で正規化    　※総電力で正規化したら間違った振幅スペクトルが得られてしまった
			for(i=0;i<N;i++){
				// Hf_multifd[i]=Hf_multifd[i]/sqrt(Psum1/N);		//H(f)は電圧なので，平均電力のルートで正規化
				Hf_multifd[i]=Hf_multifd[i]/sqrt(1.0/(double)N);		//☆☆☆　H(f)の平均電力を1Wに正規化しなければうまくいった
			}														//☆☆☆の部分をやってなくてもQPSKではなぜかうまくいった，16QAMでは爆死
			/// GI insertion ////////////////////////////////////////////////////////////////////////////
			GI_insertion(st,st_GI);
			
	
			/*マルチパス通信路通過時の考え方
			Delay=指導書で言う直接波，遅延波1，遅延波2...(i=0→直接波,i=1→遅延波1,i=2,遅延波2)
			st3=指導書で言う「受信波」→1サンプル毎に考えて（指導書の図でいうと縦で見る）ブロックを直接波,遅延波1,遅延波2...のすべてを足したものが受信波，プログラムで言う+=の部分
			ht_multifd=指導書で言う図8-1のインパルス応答，直接波の時はh0を，遅延波1の時はh1を，遅延波2の時はh2をかけて畳み込み和をとる
			*/
			///   st_GIがマルチパス通信路を通過　///////////////////////////////////////////////////////////////
			for(i=0; i<(K+GI)*s; i++)			
				st3[i]=0.0+0.0*I;		//初期化
				
			for(i=0; i<L; i++){		//直接波1個＋遅延波(L-1)個
				for(j=0; j<(K+GI)*s; j++){		
					if(j<i*s)	//i=0の時は直接波なのでこのifには入らない，遅延波1なら1サブキャリア分前回の送信信号tmpをs点コピー
						Delay[j]=tmp[j+N-(i*s)];	//遅延波i番目は前シンボルtmpのケツi*s個コピー	
					else							
						Delay[j]= st_GI[j-(i*s)];	//前シンボルの遅延をコピーし残りをst_GIをst_GI[0]から順番にコピーしていく		
						
					st3[j]+=Delay[j]*ht_multifd[i*s];	//直接波*h1+遅延波1*h2+遅延波2*h3+・・・・+遅延波(L-1)*hLとして，畳み込み和
				
			
				}			//mpfのTs毎の値を取り出し，L回分畳み込み和 i=の時は最初のインパルス，i=1の時は次のインパルス
			}
			
			
			//　次シンボル遅延波ように送信信号をコピー　//////////////////////////////////////////////////////////////
			for(i=0; i<N; i++){
				tmp[i]=st[i];		//次シンボルの遅延波に含まれるのでtmpに一時的に送信信号を格納
			}
			
			///////　ガウス雑音付加　　////////////////////////////////////////////////////////////////////////
			boxmuller((K+GI)*s,awgn,0.0,sigma);		//　AWGNの生成　
			for(i=0;i<(K+GI)*s;i++){
				rt1[i]=st3[i]+awgn[i];		//rt(t)=st(t)+AWGN(t)　
				// rt1[i]=st3[i];			//ガウス雑音を0にした場合→BERは0になる確認
			}
			
			////// GI removal　 ///////////////////////////////////////
			GI_removal(rt1,rt2);			
			
			////// FFT r2(t)→R2(f)　/////////////////////////////////////////////////////////////
			FFT(rt2,rF2);
			
			
			////// FDE(ZF criterion) Rhat(f)=R2(f)/H(f) /////////////////////////////////////////////////////////////
			for(i=0;i<N;i++){
				// rhatF[i]=rF2[i]/Hf_multifd[i];		//ZF等化
				rhatF[i]=rF2[i]*(creal(Hf_multifd[i])-I*cimag(Hf_multifd[i]))/(cabs(Hf_multifd[i])*cabs(Hf_multifd[i])+sigma/average2(st));
			}
			
			//////  P/S変換　Rhat(f)→R(f)　///////////////////////////////////////////////////////////////////////////  
			PS(rhatF,rF);
	
			
			
			//// 誤り判定(先頭点のみ)  ////////////////////////////////////////////////////////////////////////
			for(j=0; j<K; j++){
				if(creal(sF2[j*s])>=0 && creal(rF[j*s])<0) ci++;		
				if(creal(sF2[j*s])<0 && creal(rF[j*s])>=0) ci++;
				if(cimag(sF2[j*s])>=0 && cimag(rF[j*s])<0) cq++;	
				if(cimag(sF2[j*s])<0 && cimag(rF[j*s])>=0) cq++;
			}	
		}
	ber=(ci+cq)/(bit*(double)frame*K);
	printf("%d\t%e\n",EbN0,ber);
	}
	

	
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


/////////////////////////////  16QAM信号の生成　　///////////////////////////
void QAM16(double complex QAM[]){
	int i,j;
	double b1,b2,b3,b4;		//乱数用変数
	double amp;				//振幅A
	double s1I,s1Q;			//I軸，Q軸の振幅
	
	////////////// 16QAM信号生成   ////////////////////
	amp=1/sqrt(10.0*s);		//振幅Aを算出
	
	for(i=0;i<K;i++){
		b1=rand()%2;	b2=rand()%2;		//I軸の乱数
		b3=rand()%2;	b4=rand()%2;		//Q軸の乱数
		
		////////　　同相成分決定   ////////////
		if(b1==0 && b2==0) s1I=3*amp;		//(00)なら同相成分は3A
		if(b1==0 && b2==1) s1I=amp;			//(01)なら同相成分はA
		if(b1==1 && b2==1) s1I=-amp;		//(11)なら同相成分は-A
		if(b1==1 && b2==0) s1I=-3*amp;		//(10)なら同相成分は-3A
		
		////////　　直交成分決定   ////////////
		if(b3==0 && b4==0) s1Q=3*amp;		//(00)なら直交成分は3A
		if(b3==0 && b4==1) s1Q=amp;			//(01)なら直交成分はA
		if(b3==1 && b4==1) s1Q=-amp;		//(11)なら直交成分は-A
		if(b3==1 && b4==0) s1Q=-3*amp;		//(10)なら直交成分は-3A
		
		for(j=0;j<s;j++)
			QAM[i*s+j]=s1I+I*s1Q;			//1シンボルごとにOS比分コピー
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

///////////////////////GI insertion //////////////////////////////////////
void GI_insertion(double complex data[],double complex data_GI[]){
	int i,j;

	for(j=0;j<(K+GI);j++){	//サブキャリア数+GI長だけ繰り返す
		if(j<GI){		//GI長の数だけこのifの中に入る
			for(i=0;i<s;i++)
				data_GI[i+j*s]=data[s*(K-GI)+j*s+i];	//GI長×s点のデータコピー
		}else{
			for(i=0;i<s;i++)
				data_GI[i+j*s]=data[i+(j-GI)*s];		//data_GIの残りにｓｔ[N]のデータをコピー m=GIとなってるから(m-GI)とする
		}		
	}	
	
	//GI挿入正しいか確認
	/*  putchar('\n');
	printf("番号\tGI挿入後の信号data_GI\t\tGI挿入前の信号data\n");
	for(i=0;i<(K+GI)*s;i++){
		printf("%d番目\t%f\t%f\t",i,creal(data_GI[i]),cimag(data_GI[i]));
		if(i>=(GI*s)){
			printf("%f\t%f\n",creal(data[i-(GI*s)]),cimag(data[i-(GI*s)]));	
		}else{
			printf("%f\t%f\n",0.0,0.0);		//dataの足りない部分を0入れておく
		}
	}  *///これによってdata_GIの最初のGi*s個分dataの後ろのデータがコピーできていることが確認できた．
	
}

////////////////////////// GI removal　 ///////////////////////////////////////
void GI_removal(double complex data_GI[],double complex data[]){
	int i,j;
	
	for(i=0;i<K;i++){
		for(j=0;j<s;j++)
			data[j+i*s]=data_GI[j+(i+GI)*s];		//rt1の最初のGI*s分を捨てて，それ以降のN個のデータをコピーする
	}
			
}
			 
	////////////////////////////  信号の平均電力を求める関数 //////////////////////////////////////
double average2(double complex data[]){
	int i;
	double sum=0.0;
	
	for(i=0;i<N;i++)
		sum+=cabs(data[i])*cabs(data[i]);
	
	return sum/N;
}		

