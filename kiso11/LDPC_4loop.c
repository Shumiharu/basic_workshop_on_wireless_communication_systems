/*/////////////////////  タナーグラフ内の４ループの最小値と平均値算出　　///////////////////////////////

<目的>
・4ループの数が符号長Nや重み（ｊ，ｋ）でどのように変化するか調査する
・4ループの数がその条件下で最も小さい検査行列Hを書き出し保存する

*/

#include <stdio.h>
#include <stdlib.h>		//乱数使用時のプログラムsv5
#include <time.h>		//time関数(乱数の種を設定するため)使用時のプログラム
#include <math.h>		//数学的関数使用時のプログラム
#include <complex.h>	//複素数使用時のプログラム



#define s 1						//オーバーサンプル比


//LDPC符号用

#define column_weight 3		//列重みj(小さいほう)
#define row_weight 6		//行重みk(大きいほう)
#define N 512				//符号長N
#define M (N*column_weight/row_weight)		//パリティビット長M


#define lmax 40		//最大繰り返し数
#define choise 1	//1=検査行列Hを自分で作成しBERを出す，2=検査行列Hを先輩の見つけた性能の良いHをそのまま使いBERを出す

//BER特性用↓↓↓↓
#define bit 2 				//1シンボル当たりのビット数 BPSKの場合は1
#define BnTs 1   				
#define first 100000			//一旦雑にBER出す用のパラメータ


////////////////////////関数プロトタイプ宣言/////////////////////////////////
void QPSK(int nd, double complex base[]);
void inter(double complex s1[],double complex s2[]);		//先頭サンプルのみ抜き出す
void boxmuller(int num,double complex noise[],double mean,double variance);
void QPSK_bit(int nd,double complex input[],int output[]);
void bit_QPSK(int nd, int input[], double complex output[]);
int arg_min(int n1, int n2);
void matrix_prod(int n1,int n2,int n3,double complex mx1[n1][n2],double complex mx2[n2][n3],double complex mxout[n1][n3]);
void inv_matrix(int num,int nd,double complex IN[num][num][nd],double complex OUT[num][num][nd]);

/////////// 順列(permutation)0～num-1までのnum個の数字をランダムに入れ替える関数 /////////////
void shuffle(int num, int perm[]){
	int i,j;
	int tmp;
	
	//順列に0～num-1の整数を格納
	for(i=0; i<num; i++)
		perm[i]=i;
	
	//順列の要素，整数を並び替える
	for(i=0; i<num; i++){
		j=rand()%num;		//余り，0～numまでの値がjに格納される
		tmp=perm[i];		//一時的にコピー
		perm[i]=perm[j];	//要素iと要素jを入れ替える
		perm[j]=tmp;		//要素jにtmpを入れる
	}
}

////////////////　nr×nc行列matrixにおいて，指定した2つの行n1と行n2を入れ替える関数 ////////////////
void row_swap(int n1,int n2,int nr,int nc,int matrix[nr][nc]){
	int i,j,k;
	int tmp[nc];		//行列matrixの指定行n1の要素をコピー
	
	//行n1をtmpにコピー
	for(i=0; i<nc; i++){
		tmp[i]=matrix[n1][i];
	}
	
	for(i=0; i<nc; i++){
		matrix[n1][i]=matrix[n2][i];
		matrix[n2][i]=tmp[i];
	}
}

////////////////　nr×nc行列matrixにおいて，指定した2つの列n1と列n2を入れ替える関数 ////////////////
void column_swap(int n1,int n2,int nr,int nc,int matrix[nr][nc]){
	int i,j,k;
	int tmp[nr];		//行列matrixの指定列n1の要素をコピー
	
	//列n1をtmpにコピー
	for(i=0; i<nr; i++){
		tmp[i]=matrix[i][n1];
	}
	
	for(i=0; i<nr; i++){
		matrix[i][n1]=matrix[i][n2];
		matrix[i][n2]=tmp[i];
	}
}

////////////////　行数または列数numの行列において，行または列2つの排他的論理和をとる関数 ////////////////
//			num:行数または列数    data1をdata2に加える(=XOR)
void array_XOR(int num, int data1[num], int data2[num]){
	int i,j;
	
	for(i=0; i<num; i++)
		data2[i]=data2[i]^data1[i];
	
}

double Gallager_function(double x){
	double fx;
	
	if(x<0.0000001){
		fx=(log((exp(0.0000001)+1.0)/(exp(0.0000001)-1.0)));
	}else if(x>30.0){
		fx=(log((exp(30.0)+1.0)/(exp(30.0)-1.0)));
	}else{
		fx=(log((exp(x)+1.0)/(exp(x)-1.0)));
	}
	return fx;
}

double sign(double x){
	
	if(x<0)
		return -1.0;
	else
		return 1.0;	
	
}



int main(void)
{
	
	
	
///////////////////////      乱数の種を設定     ///////////////
	srand((unsigned int)time(NULL)); 
	
	int i,j,k,l,m,o;
	int loop;

	// int H[M][N];		//検査行列H(M×N)
	int (*H)[N]=(int(*)[N])malloc(M*N*sizeof(int));		//検査行列H(M×N)
	// int tmp_H[M][N];
	int (*tmp_H)[N]=(int(*)[N])malloc(M*N*sizeof(int));		//検査行列H(M×N)

////////choise=1 : 検査行列Hの生成 //////////////////////////////////////////////////////////////////////////////////////////////////////

	
	// int random_column[N];		//検査行列Hの2段目以降の要素をランダムに入れ替えるため
	int (*random_column)=(int(*))malloc(N*sizeof(int));
	int c_4loop;		//4ループの数
	int sum_4loop;
	double min_4loop;
	
	
	printf("符号長N=%d, パリティビット長M=%d, 符号化率R=%f\n",N,M,(N-M)/(double)N);
	
	sum_4loop=0;
	for(loop=0; loop<first; loop++){
		c_4loop=0;		//初期化
		
		//検査行列Hの1段目のサブブロックに階段状に1を格納
		k=0;
		for(l=0;l<N/row_weight;l++){	//N/row_weight回
			for(i=row_weight*l;i<row_weight*(l+1);i++){	//行重みkごとに同じ操作を繰り替えす
				for(j=0;j<M/column_weight;j++){		//1段目のサブブロック作成，縦向きに順に作成していく
					H[j][i]=(j==k)? 1:0;		//j=k(k段目にのみ1)，それ以外は0を格納
				}
			}
			k++;		//1を階段状に入れるためのk++
		}
		
		//検査行列Hの2段目以降のサブブロックの格納
		for(i=1; i<column_weight; i++){
			shuffle(N,random_column);		//列をランダムに入れ替えるため
			for(j=0; j<N; j++){
				for(k=0; k<M/column_weight; k++){
					H[(M/column_weight)*i+k][j]=H[k][random_column[j]];
				}//サブブロックの行数kループ　end
			}//検査行列Hの列数j end
		}//1段目以降のサブブロックiループ end 
		
		//行入れ替えのためのランダム列の確認
		// putchar('\n');
		// for(i=0; i<N; i++)
			// printf("%d ",random_column[i]);	
		
		//検査行列Hの確認
		// putchar('\n');
		// for(i=0;i<M;i++){
			// putchar('\n');
			// for(j=0;j<N;j++)
				// printf("%d ",H[i][j]);
		// }
		// putchar('\n');
		
		
		////4ループ数探査　//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		for(i=0; i<M; i++){	//行シフト
			for(j=0; j<N; j++){	//列シフト
				if(H[i][j]==1){	//成分が1ならば
					for(k=j+1; k<N; k++){	//今見ている列から右にシフト
						if(H[i][k]==1){	//今見ている行iにおいて1が列jと列kにあるなら
							for(m=i+1; m<M; m++){	//列kの行iから下に行シフト
								if(H[m][k]==1){
									for(o=k-1; o>=0; o--){	//列kから左に列シフト
										if(H[m][o]==1&& o==j){	//H[m][o]=1でかつ最初の注目した列jと一致していたら
											c_4loop++;	//4ループカウント
											// printf(" (%d,%d) (%d,%d) (%d,%d) (%d,%d) 4loop=%d \n",i,j,i,k,m,k,m,o,c_4loop);
										}
									}//o end
								}
							}//m end
						}
					}//k end
				}
			}//j end
		}//i end
			
		//4ループの最小値の初期化
		if(loop==0)	min_4loop=c_4loop;
		
		//最小値の更新とその時の検査行列のコピー
		if(c_4loop<min_4loop){
			min_4loop=c_4loop;	//最小値の更新
			for(i=0; i<M; i++){
				for(j=0; j<N; j++)
					tmp_H[i][j]=H[i][j];	//4ループ最小の検査行列Hをコピー
			}
		}				
		
		//4ループの合計値の更新
		sum_4loop+=c_4loop;
		
		
		
		
		printf("loop=%d 符号長=%d 4ループの最小値=%.0f  4ループの平均[個]=%f\r",loop,N,min_4loop,(double)sum_4loop/(loop+1));
	}//	ループn end
	
	
	///////結果の表示////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	printf("loop=%d 符号長=%d 4ループの最小値=%.0f  4ループの平均[個]=%f\n",loop,N,min_4loop,(double)sum_4loop/loop);
	
	
	///4ループ最小の検査行列の書き込み /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	FILE *fp;//ファイル作成
	char fname[100];
	
	//4ループ最小の検査行列の書き込み
	sprintf(fname, "4loopmin_Hmatrix,N=%d,M=%d,j=%d,k=%d.txt",N,M,column_weight,row_weight);
	fp=fopen(fname,"w");	//書き込みのwrite
	if(fp==NULL){
		printf("%sERROR!\n",fname);
		exit(EXIT_FAILURE);
	}
	for(i=0; i<M; i++){//検査行列Hの書き込み
		for(j=0; j<N; j++){
			fprintf(fp,"%d ", tmp_H[i][j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	printf("%s Finish Writting\n",fname);
	
	
	//メモリ開放
	free(H);			free(tmp_H);			free(random_column);
	

	
	return 0;
}


/////////////   複素ベースバンド信号(と搬送帯域信号)の生成　　　//////////////////////
void QPSK(int nd, double complex base[]){
		int i,j;
		double tmp_I;	//同相成分Iを一時的に格納
		double tmp_Q;	//直交成分Qを一時的に格納
		
		for(i=0; i<nd; i++){
			tmp_I=1/(sqrt(2*s))*(-2*(rand()%2)+1);	//電圧の大きさが1/√2s , 乱数が0なら１で1ならー1	
			tmp_Q=1/(sqrt(2*s))*(-2*(rand()%2)+1);	//電圧の大きさが1/√2s , 乱数が0なら１で1ならー1
			
			for(j=0; j<s; j++)
				base[s*i+j]=tmp_I+tmp_Q*I;				//複素ベースバンド信号u(t)=I(t)+jQ(t)
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


/////////////////// QPSK複素信号をビット系列を変換する関数 //////////////////
//             入力のデータ数     QPSK複素信号(nd)   ビット系列(nd*bit)
void QPSK_bit(int nd,double complex input[],int output[]){
	int i,j;
	
	for(i=0;i<nd;i++){
		if(creal(input[i])>0 && cimag(input[i])>0){		//第1象限(00)
			output[i*2]=0;
			output[i*2+1]=0;
		}
		if(creal(input[i])<0 && cimag(input[i])>0){		//第2象限(10)
			output[i*2]=1;
			output[i*2+1]=0;
		}
		if(creal(input[i])<0 && cimag(input[i])<0){		//第3象限(11)
			output[i*2]=1;
			output[i*2+1]=1;
		}
		if(creal(input[i])>0 && cimag(input[i])<0){		//第4象限(01)
			output[i*2]=0;
			output[i*2+1]=1;
		}	
	}
}

/////////////////////　ビット列をQPSK信号に変換する関数  //////////////////////
//			    ビット長   ビット系列(nd)        QPSK信号(nd/2) 
void bit_QPSK(int nd, int input[], double complex output[]){
	int i,j;
	
	for(i=0;i<nd/2; i++){
		if(input[i*2]==0 && input[i*2+1]==0)
			output[i]=1/(sqrt(2*s))+I*1/(sqrt(2*s));	//第1象限(00)=A+jA
		
		if(input[i*2]==1 && input[i*2+1]==0)
			output[i]=-1/(sqrt(2*s))+I*1/(sqrt(2*s));	//第2象限(10)=-A+jA
		
		if(input[i*2]==1 && input[i*2+1]==1)
			output[i]=-1/(sqrt(2*s))-I*1/(sqrt(2*s));	//第1象限(11)=-A-jA
		
		if(input[i*2]==0 && input[i*2+1]==1)
			output[i]=1/(sqrt(2*s))-I*1/(sqrt(2*s));	//第1象限(01)=A-jA
	}
}


int arg_min(int n1, int n2){
	if(n1<=n2)
		return n1;		
	else
		return n2;
}
	

/////////////////  num1×num2とnum2×num3行列の掛け算 出力はnum1×num3の行列 ///////////////////////////////// 
void matrix_prod(int n1,int n2,int n3,double complex mx1[n1][n2],double complex mx2[n2][n3],double complex mxout[n1][n3]){
	int i,j,k,l;
	
	
	for(j=0;j<n1;j++){
		for(k=0;k<n3;k++){
			mxout[j][k]=0.0+0.0*I;
		}
	}

	
	for(j=0;j<n1;j++){
		for(k=0;k<n3;k++){
			for(l=0;l<n2;l++){
				mxout[j][k]+=mx1[j][l]*mx2[l][k];
			}
		}
	}

}

//////////////////////////////////// num×numの行列の逆行列を求める関数 ////////////////////////////
//                       nd:データ数
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














