/*
 * rand_C.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 */

/*	rand_C.cpp

1. 简介： 
	提供随机数生成及随机树结构初始化的辅助功能，支持 MCMC 过程中的随机抽样。

2.背景说明： 
	该文件为整个算法提供随机化支持。在 MCMC 迭代中，它用于决定是否接受新的参数提案（Moves）、
	选择变动类型以及初始化随机的突变树或二叉样本树，对应于随机搜索算法的初始化和状态转移采样步骤。

3.主要内容：
	initRand()：初始化随机数生成器种子 
	getRandTreeCode()：生成长度为 n-1 的随机 Prüfer 序列，用于构建随机有根树 
	changeBeta()：根据给定概率决定当前 MCMC 步是否更新错误率参数 
	sampleRandomMove()：根据移动概率分布随机挑选一种树结构变动类型 
	samplingByProb()：通用概率抽样函数，返回是否命中给定概率 
	sampleTwoElementsWithoutReplacement()：从 n 个元素中无放回地随机抽取两个不同的元素 
	pickRandomNumber()：从 [0, n-1] 区间内随机抽取一个整数 sample_0_1()：生成 [0, 1] 之间的随机浮点数 
	getElemFromQueue()：从向量中获取指定索引的元素并将其与末尾交换以优化删除 
	getRandomBinaryTree()：通过随机合并过程构建一棵包含 m 个叶节点的随机二叉树父节点向量

4.输入：
	prob：变动概率或接受概率值（double 类型或 vector 容器） 
	n：突变数量或节点上限 
	m：叶节点（细胞）数量

5.输出：
	code/result/leafsAndInnerNodesParents：动态分配的整型数组指针，包含随机生成的树编码、索引对或父节点向量 
	true/false：布尔值，指示抽样或决策结果

6.注意事项：
	算法假设：
		随机性基于标准库的 rand() 函数，假设调用者已通过 initRand() 完成初始化；构建二叉树时假设内节点索引从 m 开始。 

	使用时的重要限制或潜在问题：
		使用简单的模运算 (rand() % n) 在 n 较大时可能存在轻微的均匀性偏差；返回的指针指向动态分配的内存，需由调用者负责释放（delete []）。
*/

//#include <string>
//#include <random>
#include "rand.h"
#include <iostream>
#include <random>
#include <stdlib.h>
#include <time.h>
#include "matrices.h"

using namespace std;


/*****    functions for sampling random numbers inside C++  *****/
void initRand(){
	time_t t;
	time(&t);
	srand((unsigned int)t);              // initialize random number generator
	//srand(1);
}


/* This function gets a number of nodes n, and creates a random pruefer code for a rooted tree with n+1 nodes (root is always node n+1) */
int* getRandTreeCode(int n){                // as usual n is the number of mutations

	int nodes = n+1;                        // #nodes = n mutations plus root (wildtype)
	int codeLength = nodes-2;
	int* code = new int[codeLength];
	for(int i=0; i<codeLength; i++){
		code[i] = rand() % nodes;
	}
	return code;
}

bool changeBeta(double prob){
	 double percent = (rand() % 100)+1;    // between 1 and 100
	 if(percent <= prob*100){
		 return true;
	 }
	 return false;
}

int sampleRandomMove(std::vector<double> prob){ // picks randomly one of the tree moves based on the move probabilities

    double percent = (rand() % 100)+1;    // between 1 and 100
    double probSum = prob[1];
    for(int i=1; i<prob.size()-1; i++){    // start at index 1; the probability at prob[0] is for changing the error rate (which is treated separately)
        if(percent <= probSum*100){
          return i;
        }
        probSum += prob[i+1];
    }
    return prob.size()-1;
}


bool samplingByProb(double prob){
	double percent = rand() % 100;
	if(percent <= prob*100){
		return true;
	}
	return false;
}


int* sampleTwoElementsWithoutReplacement(int n){

    int* result = new int[2];
	  result[0] = rand() % n;
	  result[1] = result[0];
    while(result[0]==result[1]){
      result[1] = rand() % n;
    }
	  return result;
}

int pickRandomNumber(int n){

    return (rand() % n);
}

double sample_0_1(){

  //return (((double) rand()+0.5) / ((RAND_MAX+1)));
  return ((double) rand() / RAND_MAX);
}

int getElemFromQueue(int index, std::vector<int> queue){
	int elem = queue.at(index);
	if (index != queue.size() - 1)
	{
		queue[index] = std::move(queue.back());
	}

	//cout << queue.size() << " elements in queue in subroutine\n";
	return elem;
}

// This creates the parent vector of a random binary tree. Entries 0...m-1 are for the leafs.
// Entries m....2m-3 are for the inner nodes except the root, the root has index 2m-2 which has no parent
// and therefore has no entry in the parent vector
int* getRandomBinaryTree(int m){
	int parentCount = (2*m)-2;     // the m leafs have a parent and so have m-2 of the m-1 inner nodes
	int* leafsAndInnerNodesParents = init_intArray(parentCount, -1);

	std::vector<int> queue;
	for(int i=0; i<m; i++){queue.push_back(i);}   // add the m leafs to the queue
	//cout << queue.size() << " elements in queue\n";
	int innerNodeCount = m;
	while(queue.size()>1){
		int pos = pickRandomNumber(queue.size());
		int child1 = queue.at(pos);
		if (pos != queue.size() - 1){queue[pos] = std::move(queue.back());}
		queue.pop_back();
		//cout << queue.size() << " elements in queue\n";

		pos = pickRandomNumber(queue.size());
		int child2 = queue.at(pos);
		if (pos != queue.size() - 1){queue[pos] = std::move(queue.back());}
		queue.pop_back();
		//cout << queue.size() << " elements in queue\n";

		leafsAndInnerNodesParents[child1] = innerNodeCount;
		leafsAndInnerNodesParents[child2] = innerNodeCount;
		queue.push_back(innerNodeCount);
		innerNodeCount++;
	}
	return leafsAndInnerNodesParents;
}

