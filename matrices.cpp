/*
 * matrices.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 */

 /*	matrices.cpp

1. 简介： 
	提供矩阵与数组的基础操作函数库，支持生物信息学算法中的数据结构处理。

2. 背景说明： 
	该文件为 SCITE 项目的基础工具模块，负责处理肿瘤进化推断中涉及的突变矩阵、祖先矩阵（Ancestry Matrix）以及树结构转换。
	特别是 ancMatrixToParVector 函数，在将推断的树矩阵转换为父节点向量表示法（MCMC 搜索的核心步骤）中起到关键作用。

3. 主要内容： 
	getMaxEntry()：返回双精度浮点数组中的最大值 
	sumMatrices()：将两个整数矩阵相加并返回结果 
	transposeMatrix()：对整数矩阵进行转置操作 
	addToMatrix()：将第二个整数矩阵累加到第一个矩阵上 
	ancMatrixToParVector()：将祖先矩阵转换为树的父节点向量表示 
	allocate_doubleMatrix/intMatrix/boolMatrix()：动态分配不同类型的二维数组内存 
	init_intArray/doubleArray/boolArray()：初始化一维数组并填充指定初值 
	init_doubleMatrix/intMatrix/boolMatrix()：分配并初始化二维矩阵 
	reset_intMatrix()：将现有整数矩阵的所有条目重置为指定值 
	delete_3D_intMatrix()：释放三维整数矩阵的内存 
	free_boolMatrix/intMatrix/doubleMatrix()：释放二维矩阵的内存 
	deepCopy_boolMatrix/intMatrix/doubleMatrix()：执行二维矩阵的深拷贝 
	deepCopy_intArray/doubleArray()：执行一维数组的深拷贝 
	identical_boolMatrices()：比较两个布尔矩阵是否完全一致 
	print_boolMatrix/doubleMatrix/intMatrix()：格式化打印矩阵内容 
	print_intArray()：格式化打印整数数组内容

4. 输入： 
	matrix/first/second：待处理的输入矩阵指针 
	n：矩阵的行数或数组长度 
	m：矩阵的列数 
	value：初始化时使用的填充值 
	anc：表示节点间祖先关系的布尔矩阵

5. 输出： 
	transposed：转置后的新矩阵指针 
	parVec：转换得到的父节点向量数组指针 
	deepCopy：拷贝生成的新数组/矩阵指针 
	maxEntry：查找到的最大数值

6. 注意事项： 
	算法假设：
		ancMatrixToParVector 函数假设输入的祖先矩阵能够代表一个合法的树结构，其中节点编号 0 到 n-1 为突变，n 为根节点。 
	使用时的重要限制或潜在问题：
		内存管理采用 C++ 标准的 new[] 和 delete[]，用户需确保调用对应的 free/delete 函数以避免内存泄漏；矩阵分配采用了连续内
		存块优化（matrix[i] = matrix[i-1] + m），释放时需严格遵守配套的 free 函数逻辑。
 */
#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include "matrices.h"

using namespace std;

/*****    basic functions on 1D and 2D arrays   *****/

double getMaxEntry(double* array, int n){
	double maxEntry = -DBL_MAX;
	for(int i=0; i<n; i++){
		maxEntry = max(maxEntry, array[i]);
	}
	return maxEntry;
}

int** sumMatrices(int** first, int** second, int n, int m){
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			first[i][j] += second[i][j];
		}
	}
	return first;
}

int ** transposeMatrix(int** matrix, int n, int m){
	int ** transposed = allocate_intMatrix(m, n);
	for(int i=0; i<m; i++){
		for(int j=0; j<n;j++){
			transposed[i][j] = matrix[j][i];
		}
	}
	return transposed;
}

void addToMatrix(int** first, int** second, int n, int m){
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			first[i][j] += second[i][j];
		}
	}
}

int* ancMatrixToParVector(bool** anc, int n){
	int* parVec = new int[n];
	for(int i=0; i<n; i++){
		parVec[i] = n;
	}
	for(int i=0; i<n; i++){
		for(int k=0; k<n; k++){
			if(k!=i && anc[k][i]==true){  // k is true ancestor of i
				bool cand = true;
				for(int l=0; l<n; l++){
					if(l!=i && l!=k && anc[l][i]==true && anc[k][l]==true){   // k is ancestor of l, and l is ancestor of i
						cand = false;                                        // k is no longer candidate for being parent of i
						break;
					}
				}
				if(cand==true){           // no ancestor of i is descendant of k -> k is parent of i
						parVec[i] = k;
				}
			}

		}
	}
	return parVec;
}

/*   allocation  */
double** allocate_doubleMatrix(int n, int m){

    double** matrix = new double*[n];
    matrix[0] = new double[n*m];
	  for (int i=1; i<n; ++i)
    {
        matrix[i] = matrix[i-1] + m;
    }
    return matrix;
}

int** allocate_intMatrix(int n, int m){

    int** matrix = new int*[n];
    matrix[0] = new int[n*m];
    for (int i=1; i<n; ++i)
    {
        matrix[i] = matrix[i-1] + m;
    }
    return matrix;
}

bool** allocate_boolMatrix(int n, int m){

    bool** matrix = new bool*[n];
    matrix[0] = new bool[n*m];
    for (int i=1; i<n; ++i)
    {
        matrix[i] = matrix[i-1] + m;
    }
    return matrix;
}


/*  initialization  */

int* init_intArray(int n, int value){
	int* array = new int[n];
	for(int i=0; i<n; i++){
		array[i] = value;
	}
	return array;
}

double* init_doubleArray(int n, double value){
	double* array = new double[n];
	for(int i=0; i<n; i++){
		array[i] = value;
	}
	return array;
}

bool* init_boolArray(int n, bool value){
	bool* array = new bool[n];
	for(int i=0; i<n; i++){
		array[i] = value;
	}
	return array;
}

double** init_doubleMatrix(int n, int m, double value){

	  double** matrix = allocate_doubleMatrix(n, m);     // allocate

    for (int i=0; i<n; ++i)             // initialize
    {
         for (int j=0; j<m; ++j)
      {
        	matrix[i][j] = value;
    	}
    }
    return matrix;
}

int** init_intMatrix(int n, int m, int value){

    int** matrix = allocate_intMatrix(n, m);  // allocate

    for (int i=0; i<n; ++i)            // initialize
    {
         for (int j=0; j<m; ++j)
      {
        	matrix[i][j] = value;
    	}
    }
    return matrix;
}

void reset_intMatrix(int** matrix, int n, int m, int value){

    for (int i=0; i<n; ++i)            // initialize
    {
         for (int j=0; j<m; ++j)
      {
        	matrix[i][j] = value;
    	}
    }
}


bool** init_boolMatrix(int n, int m, bool value){

    bool** matrix = allocate_boolMatrix(n, m);     // allocate

    for (int i=0; i<n; ++i)             // initialize
    {
         for (int j=0; j<m; ++j)
      {
        	matrix[i][j] = value;
    	}
    }
    return matrix;
}


/*  deallocation  */

void delete_3D_intMatrix(int*** matrix, int n){
	for(int i=0; i<n; i++){
			delete [] matrix[i][0];
			delete [] matrix[i];
		}
		delete [] matrix;
}


void free_boolMatrix(bool** matrix){
    delete [] matrix[0];
    delete [] matrix;
}

void free_intMatrix(int** matrix){
    delete [] matrix[0];
    delete [] matrix;
}

void free_doubleMatrix(double** matrix){
    delete [] matrix[0];
    delete [] matrix;
}


/*  deep copying  */

bool** deepCopy_boolMatrix(bool** matrix, int n, int m){
    bool** deepCopy = init_boolMatrix(n,m, false);
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
    	      deepCopy[i][j] = matrix[i][j];
	      }
    }
    return deepCopy;
}

int** deepCopy_intMatrix(int** matrix, int n, int m){
    int** deepCopy = init_intMatrix(n,m, -1);
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
    	      deepCopy[i][j] = matrix[i][j];
	      }
    }
    return deepCopy;
}

double** deepCopy_doubleMatrix(double** matrix, int n, int m){
    double** deepCopy = init_doubleMatrix(n,m, -1);
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
    	      deepCopy[i][j] = matrix[i][j];
	      }
    }
    return deepCopy;
}

int* deepCopy_intArray(int* array, int n){
	  int* deepCopy = new int[n];
	  for (int i=0; i<n; ++i)
    {
        deepCopy[i] = array[i];
    }
    return deepCopy;
}

double* deepCopy_doubleArray(double* array, int n){
	double* deepCopy = new double[n];
	for (int i=0; i<n; ++i)
    {
        deepCopy[i] = array[i];
    }
    return deepCopy;
}

bool identical_boolMatrices(bool** first, bool** second, int n, int m){
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			//cout << "[" << i << "," << j << "] ";
			if(first[i][j] != second[i][j]){
				cout << "matrices differ!!!!!!!!!!!!!!!!\n";
				getchar();
				return false;
			}
		}
		//cout << "\n";
	}
	return true;
}

/*  printing  */

void print_boolMatrix(bool** array, int n, int m){
	  for(int i=0; i<n; i++){
		    for(int j=0; j<m; j++){
			      cout << array[i][j] << " ";
		    }
		    cout << "\n";
	  }
}

void print_doubleMatrix(double** matrix, int n, int m){
	  for(int i=0; i<n; i++){
  		  for(int j=0; j<m; j++){
  			    cout << matrix[i][j] << "\t";
  		  }
  		  cout << "\n";
  	}
}


void print_intArray(int* array, int n){
	for (int i=0; i<n; ++i)
    {
       cout << array[i] << " ";
    }
    cout << "\n";
}

void print_intMatrix(int** matrix, int n, int m, char del){
	  for(int i=0; i<n; i++){
  		  for(int j=0; j<m; j++){
  			    cout << matrix[i][j] <<  del;
  		  }
  		  cout << "\n";
  	}
}



