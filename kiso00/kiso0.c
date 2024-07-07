#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define T_s 0.001 /* T_s秒ごとにデータを生成する */

/* 関数プロトタイプ宣言 */
void generate_data(int symbol, int data[]); 

int main(void) {
    int i; /* for文に用いる変数 */
    srand((unsigned int)time(NULL)); /* ランダム関数の初期化 */

    int symbol; /* シンボル数（データ長）*/
    printf("シンボル数（データ長）を入力: ");
    scanf("%d", &symbol); /* シンボル数（データ数）の入力を要求 */

    int data[symbol]; /* データ */
    generate_data(symbol, data); /* 0と1のデータをランダム生成 */

    double t = 0.0; /* 経過時間 */
    printf("time[s] data\n");
    for(i = 0; i < symbol; i++) { 
        printf("%f ", t); /* 経過時間を出力 */
        printf("%d \n", data[i]); /* 生成したデータを出力 */
        t += T_s; /* 経過時間の加算 */
    }

    return 0; /* プログラムの終了 */
}

/* 0と1のデータをランダム生成 */
void generate_data(int symbol, int data[]) {
    int i;
    for(i = 0; i < symbol; i++) {
        data[i] = rand()%2; /* ランダム生成した数を2で割った余りから生成 */
    }
}
