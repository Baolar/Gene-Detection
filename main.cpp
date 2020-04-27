/*
 * 基因检测二进制问题
 * Author:Wenshuo Chen, WUST
 * mail: chenwenshuo@wust.edu.cn
 * 2020.4.27
 */

#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>
using namespace std;
const double eps = 1e-5; //精度
const double rowsum = 1;//所有为1位的强度和,调高可以解决仪器精度问题
int n, m;
double tb[500][20]; //编号-试管-权值表，tb[i][j]表示编号为i的二进制从高位到低位的第j位的权值
double csq[20];//测试试管结果
vector<double> sum; //dfs当前和
vector<vector<int>> ans; //存答案
vector<int> used; //dfs当前所用编号

/*
 * 暴力搜索枚举
 * 枚举每一个编号是否被用过
 */
void dfs(int x, vector<double> sum, vector<int> used) { //sum维护当前所选编号权值和
	bool flag = 1;
	for (int i = 0; i < sum.size(); i++) {
		if (sum[i] > csq[i + 1] + eps) return; //若存在所选编号i试管的权值和大于该编号的结果，剪枝
		else if (fabs(sum[i] - csq[i + 1]) > eps) {	///如果存在当前和和结果不在误差范围内，置flag = 0		
			flag = 0;
			break;
		}
	}
	if (flag) { //flag = 1代表当前used所选编号的试管强度和为实验结果值，将它存入答案
		ans.push_back(used);
		return;
	}
	
	if (x > n) return; //遍历完成，结束递归
	
	dfs(x + 1, sum, used); //编号x的基因不选（不突变），used不变
	for (int i = 0; i < sum.size(); i++)
		sum[i] += tb[x][i + 1];
	used.push_back(x); 
	dfs(x + 1, sum, used); //若编号x的基因突变，将编号x存入used并将编号x的二进制每一位权值加到sum上
}

/*
 * 获取二进制位为1的位数
 * 比如 5 为 101 ，返回2
 */
int get1bcnt(int x) { 
	int res = 0;
	while (x) {
		if (x & 1) res++;//最后一位与1取与，为true即最后一位是1，res++
		x >>= 1; //右移一位
	}

	return res;
}

void create_tb(int n, int m) {
	for (int i = 1; i <= n; i++) {
		int row1cnt = get1bcnt(i); //获取编号i的二进制位为1的位数

		int x = i;
		for (int j = m; j >= 1; j--) {
			if (x & 1) tb[i][j] = rowsum / row1cnt; //每一位位为1的权值就是预置和/位数
			else tb[i][j] = 0;
			x >>= 1;
		}
	}
}

void display_tb(int n, int m) {
	printf("gene\\tb");
	for (int i = 1; i <= m; i++)
		printf("%8d", i);
	printf("\n");
	for (int i = 1; i <= n; i++) {
		printf("%8d", i);
		for (int j = 1; j <= m; j++) {
			printf("%8.3f", tb[i][j]);
		}
		printf("\n");
	}
}

int main() {
	cout << "                    基因喝试管问题" << endl;
	cout << "Author: Wenshuo Chen, Wuhan University of Science and Technology" << endl;
	cout << "                      2020.4.27 " << endl << endl;
	cout << "注：请注意输入的格式否则会崩溃，测试结果的值不易太大" << endl;
	cout << "    即变异基因不要太多，最好40以内，每一个试管检测结果信号强度应该合理，不能过高" << endl << endl;
	while (1) {
		cout << "请输入待检测基因的个数:  ";
		do {
			cin >> n;
			if (n >= 500) cout << "程序最多仅仅设置了500个编号的检测基因,请重新输入:" << endl;
		} while (n >= 500);
	
		int nbak = n;
		m = 0;

		while (nbak) {
			nbak >>= 1;
			m++;
		}
		
		create_tb(n, m);
		
		cout << "最少需要试管数量\t" << m << endl << endl;
		
		cout << "创建编号i所对应试管试剂量的表格如下" << endl;
		display_tb(n, m);
		cout << "=========================================" << endl;
		cout << "请输入 " << m << " 个试管的检测结果(小数并用空格隔开):" << endl;
		for (int i = 1; i <= m; i++) cin >> csq[i];

		sum.resize(m);
		ans.clear();
		used.clear();

		dfs(1, sum,used);

		cout << "精度为 " << eps << endl;
		cout << "以下编号的基因可能产生了突变:" << endl;
		for (int i = 0; i < ans.size(); i++) {
			for (int j = 0; j < ans[i].size(); j++)
				cout << ans[i][j] << " ";
			cout << endl;
		}
	}

	return 0;
}
