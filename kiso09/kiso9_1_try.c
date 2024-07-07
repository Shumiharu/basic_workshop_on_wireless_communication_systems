#include <stdio.h>
#include <stdlib.h>		//乱数使用時のプログラム
#include <time.h>		//time関数(乱数の種を設定するため)使用時のプログラム
#include <math.h>		//数学的関数使用時のプログラム
#include <complex.h>	//複素数使用時のプログラム


#define a 1			 			//波数
#define s 1						//オーバーサンプル比
#define K 64					//サブキャリア数
#define N (s*K)					//全サンプル数
#define sigmaN 0.5				//ガウス雑音の分散
#define Ts 0.001				//1サンプル時間(OFDMでは)
#define fc a/Ts 				//搬送波周波数fc
#define T_sample (Ts/s) 		//1サンプル当たりの時間
#define delta_f (1/(Ts*K))		//Δf=1/T Δf:スペクトルの幅  T(=Ts*K):時間長

//OFDM特有のパラメータ↓↓↓↓
#define L 1						//L波マルチパス
#define dB 1					//何dB減衰
#define GI 4					//ガードインターバルのサンプル数

//MIMO特有のパラメータ
#define Nt 2					//送信アンテナNt
#define Nr 2					//受信アンテナNr


//BER特性用↓↓↓↓
#define bit 2.0 				//1シンボル当たりのビット数 QPSKの場合は2
#define BnTs 1   				
#define frame pow(10,5)			//フレーム数

//n×nの逆行列
/* 
void inverse_matrix(int num,double complex matrix1[][],double complex IN[Nt][num][N]){
	int i,j,k,l;
	double complex temp1,temp2,temp[Nt][num][N];
	//tempにINを格納
	for(l=0;l<N;l++){
		for(i=0;i<Nt;i++){
			for(j=0;j<num;j++){
				temp[i][j][l]=IN[i][j][l];
			}
		}
	}
	
	//単位行列生成
	for(l=0;l<N;l++){
		for(i=0;i<Nt;i++){
			for(j=0;j<num;j++){
				if(i==j){
					OUT[i][j][l]=1.0;		
				}
				else{
					OUT[i][j][l]=0.0;
				}
			}
		}
	}
	
	//掃き出し
	for(l=0;l<N;l++){
		for(i=0;i<num;i++){
			temp1=temp[i][i][l];//i行目でほかの行を割る
			for(j=0;j<Nt;j++){
				temp[i][j][l]=temp[i][j][l]/temp1;//単位行列を目指す
				OUT[i][j][l]=OUT[i][j][l]/temp1;//逆行列を目指す
			}
			for(k=0;k<Nt;k++){
				if(i!=k){//左側の行列が単位行列でないとき
					temp2=temp[k][i][l];//ななめの部分の数をひいて0にする
					for(j=0;j<num;j++){
						temp[k][j][l]=temp[k][j][l]-temp[i][j][l]*temp2;
						OUT[k][j][l]=OUT[k][j][l]-OUT[i][j][l]*temp2;
					}
				}
			}
		}
	}
}
 */
 
//////////////////////////////////// num×numの行列の逆行列を求める関数 ////////////////////////////
//n×nの逆行列
void inv_matrix(int num,int nd,double complex IN[num][num][nd],double complex OUT[num][num][nd]){
	int i,j,k,l;
	double complex temp1,temp2,temp[num][num][nd];
	
	//tempにINを格納
	for(l=0;l<nd;l++){
		for(i=0;i<num;i++){
			for(j=0;j<num;j++){
				temp[i][j][l]=IN[i][j][l];
			}
		}
	}
	
	//単位行列生成
	for(l=0;l<nd;l++){
		for(i=0;i<num;i++){
			for(j=0;j<num;j++){
				if(i==j){
					OUT[i][j][l]=1.0;		
				}
				else{
					OUT[i][j][l]=0.0;
				}
			}
		}
	}
	
	//掃き出し
	for(l=0;l<nd;l++){
		for(i=0;i<num;i++){
			temp1=temp[i][i][l];//i行目でほかの行を割る
			for(j=0;j<num;j++){
				temp[i][j][l]=temp[i][j][l]/temp1;//単位行列を目指す
				OUT[i][j][l]=OUT[i][j][l]/temp1;//逆行列を目指す
			}
			for(k=0;k<num;k++){
				if(i!=k){//左側の行列が単位行列でないとき
					temp2=temp[k][i][l];//ななめの部分の数をひいて0にする
					for(j=0;j<num;j++){
						temp[k][j][l]=temp[k][j][l]-temp[i][j][l]*temp2;
						OUT[k][j][l]=OUT[k][j][l]-OUT[i][j][l]*temp2;
					}
				}
			}
		}
	}
}

/////////////////////////////////  num1×num2とnum2×num3行列の掛け算  ////////////////////////////////////
/* void matrix_prod(int num1,int num2,int num3,int nd,double complex mx1[num1][num2][data_length],double complex mx2[num2][num3][data_length],double complex mxout[num1][num3][data_length]){
	int i,j,k,l;
	
	
	for(i=0;i<num1;i++){
		for(j=0;j<num3;j++)
			mxout[i][j]=0+0*I;
	}
	
	// for(j=0;j<num2;j++){
		// for(k=0;k<num2;k++){
			// for(l=0;l<num1;l++)
				// mxout[j][k]+=mx1[j][l]*mx2[l][k];
		// }
	// }
	
	
	for(j=0;j<num1;j++){
		for(k=0;k<num3;k++){
			for(l=0;l<num2;l++)
				mxout[j][k]+=mx1[j][l]*mx2[l][k];
		}
	}
	//逆行列を出力
	 putchar('\n');
	for(i=0;i<num1;i++){
		for(j=0;j<num2;j++)
			printf("%f\t%f\t%f\n",creal(mx1[i][j]),creal(mx2[i][j]),creal(mxout[i][j]));
 
	}
	

}
*/
 
//// 
//									   	データ数
void matrix_prod(int n1,int n2,int n3,int nd,double complex mx1[n1][n2][nd],double complex mx2[n2][n3][nd],double complex mxout[n1][n3][nd]){
	int i,j,k,l;
	
	for(i=0;i<nd;i++){
		for(j=0;j<n1;j++){
			for(k=0;k<n3;k++){
				mxout[j][k][i]=0.0+0.0*I;
			}
		}
	}
	
	
	
	for(i=0;i<nd;i++){
		for(j=0;j<n1;j++){
			for(k=0;k<n3;k++){
				for(l=0;l<n2;l++){
					mxout[j][k][i]+=mx1[j][l][i]*mx2[l][k][i];
				}
			}
		}
	}
}

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

int main(void)
{	
	int i,j,m; 
	
///////////////////////      乱数の種を設定     ///////////////
	srand((unsigned int)time(NULL)); 

////////////////////// 時間と周波数の設定　　//////////////////////
	double t[N]; 						//時間t    0～K*Tsまでの時間
	double f[N];						//周波数f  -Ts/s～Ts/sまでの周波数
	
	for(i=0;i<N;i++){
		t[i]=T_sample*i;
		f[i]=delta_f*(i-N/2);
	}
	
	
	
	
	double complex a7[2][2][3];
	double complex a8[2][2][2];
	
	// for(i=0;i<2;i++){
		// for(j=0;j<2;j++){
			// for(m=0;m<3;m++)
				// a7[i][j][m]=1+i+m+j;
		// }
	// }
	
	
		// inv_matrix(2,3,a7,a8);
	
	// printf("逆行列の確認\n");
	// for(i=0;i<2;i++){
		// for(j=0;j<2;j++){
			// for(m=0;m<3;m++)
			// printf("a8=%f\n",creal(a8[i][j][m]));
		// }
	// }
	
	
	//ぎょうれつの逆行列を求める関数が正しいかの確認//
	double complex a1[4][4][1]={{1+I,2,0,-1},{-1+4*I,1,2,0},{2,0,1-8*I,1},{1,-2,-1,1}};
	double complex a2[4][4][1];
	double complex a3[4][4][1];
	
	inv_matrix(4,1,a1,a2);
	matrix_prod(4,4,4,1,a1,a2,a3);
	putchar('\n');
	for(i=0;i<4;i++){
		for(j=0;j<4;j++)
			printf("a2=%f\ta2=%f\ta3=%f\ta3=%f\n",creal(a2[i][j][0]),cimag(a2[i][j][0]),creal(a3[i][j][0]),cimag(a3[i][j][0]));
	}
		
	
	//行列の積を求める関数が正しいかの確認//
	double complex a4[3][4][1]={{1,2,3,4},{2,3,4,5},{5,6,4,7}};
	double complex a5[4][1][1]={1,0,1,0};
	double complex a6[3][1][1];
	
	matrix_prod(3,4,1,1,a4,a5,a6);
	
	printf("3×4行列a4と4×1行列a5の積a6(3×1)\n");
	for(i=0;i<3;i++){
		for(j=0;j<1;j++)
			printf("a6=%f\n",creal(a6[i][j][0]));
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
			 
			






















