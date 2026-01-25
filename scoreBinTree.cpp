/*
 * scoreBinTree.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: jahnka
 */

 /*	scoreBinTree.cpp
 1. 简介：
 	计算单细胞突变数据在给定二叉样本树拓扑结构下的对数似然得分。

 2. 背景说明：
 	该文件是 SCITE 算法在处理二叉树推断时的核心评分组件。它实现了在样本演化树（Cell Tree）上寻找突变最佳放置位置的逻辑，
	对应于论文中关于计算给定树拓扑结构概率的步骤，通过动态规划或高效的树遍历来优化计算效率。

 3. 主要内容：
 	getBinSubtreeScore()：使用后序遍历方式递归计算特定突变在各个子树中（全为存在或全为缺失状态）的累积似然得分。
	getBinTreeMutScore()：通过组合子树得分，遍历树中所有可能的节点作为突变发生点，寻找单个突变的最优放置得分。
	getBinTreeScore()：计算整个突变矩阵在给定二叉树下的总最大似然得分，即所有突变最优得分之和。

 4. 输入：
 	obsMutProfiles：观测到的单细胞突变谱矩阵（$m \times n$）。
	n：突变的总数。
	m：单细胞样本的总数。
	logScores：预计算的错误对数似然矩阵。
	parent：二叉树的父节点向量表示。

 5.输出：
 	sumScore：整个模型在当前树结构下的对数似然总分（double 类型）。

 6.注意事项：
 	算法假设：
		输入必须是严格的二叉树结构；节点总数固定为 $2m-1$（$m$ 个叶节点和 $m-1$ 个内部节点）。
	使用时的重要限制或潜在问题：
		计算依赖于广度优先遍历（BFT）序来模拟后序遍历过程；突变被假设只能在树的一个位置发生（单次起源假设）。
*/

#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <iostream>
#include "limits.h"
#include "scoreBinTree.h"
#include "matrices.h"
#include "trees.h"

using namespace std;

/* computes the log likelihood for a single mutation for all subtrees of the binary tree, where the expected */
/* state of the mutation can be either absent or present in the whole subtree (passed as 'state' to the function) */
double* getBinSubtreeScore(bool state, int* bft, vector<vector<int> > &childLists, int mut, int nodeCount, int m, int** obsMutProfiles, double ** logScores){
	double* score = init_doubleArray(nodeCount, 0.0);
	for(int i=nodeCount-1; i>=0; i--){
		int node = bft[i];

		if(node < m){
			score[node] = logScores[obsMutProfiles[node][mut]][state];   // for leafs the score is just P(Dij|Eij)
		}
		else{                                                          // for inner nodes the score is the sum of the scores of the children
			if(childLists.at(node).size()!=2){
				cerr << "Error node " << node << " has " << childLists.at(node).size() << " children\n";  // tree should be binary, but isn't
			}
			score[node] = score[childLists.at(node).at(0)] + score[childLists.at(node).at(1)];
		}
	}
	return score;
}


/* Computes the best log likelihood for placing a single mutation in a given sample tree */
/* Iterates through all nodes as possible placements of the mutation to find the best one */
/* All samples below the placement of the mutation should have it, mutation can also be placed at a leaf, i.e. uniquely to the sample   */
double getBinTreeMutScore(int* bft, vector<vector<int> > &childLists, int mut, int nodeCount, int m, int** obsMutProfiles, double ** logScores){

	double bestScore = -DBL_MAX;
	double* absentScore = getBinSubtreeScore(0, bft, childLists, mut, nodeCount, m, obsMutProfiles, logScores);
	double* presentScore = getBinSubtreeScore(1, bft, childLists, mut, nodeCount, m, obsMutProfiles, logScores);

	for(int p=0; p<nodeCount; p++){
		double score = absentScore[nodeCount-1] - absentScore[p] + presentScore[p];
		bestScore = max(bestScore, score);
	}

	delete [] absentScore;
	delete [] presentScore;
	return bestScore;
}

/* Computes the maximum log likelihood of a binary tree for a given mutation matrix.  */
/* Note: No extra root necessary for binary trees */
double getBinTreeScore(int** obsMutProfiles, int n, int m, double ** logScores, int* parent){

	int nodeCount = (2*m)-1;   // number of nodes in binary tree: m leafs, m-1 inner nodes (includes already the root)
	double sumScore = 0;       // sum of maximal scores of all samples
	vector<vector<int> > childLists = getChildListFromParentVector(parent, nodeCount-1);
	int* bft = getBreadthFirstTraversal(parent, nodeCount-1);

	for(int mut=0; mut<n; mut++){                                                // sum over the optimal scores of each sample
		double score = getBinTreeMutScore(bft, childLists, mut, nodeCount, m, obsMutProfiles, logScores);
		sumScore += score;
	}

	delete [] bft;
	for(int i=0; i<childLists.size(); i++){
		childLists[i].clear();
	}
	childLists.clear();
	//cout << "score: " << sumScore << "\n";
	return sumScore;
}


//
///* Score contribution by a specific mutation when placed at a specific leaf */
///* This is the same for all trees and can be precomputed */
//double binTreeLeafScore(int** obsMutProfiles, int leaf, int mut, int m, double ** logScores){
//
//	double score = logScores[obsMutProfiles[leaf][mut]][1];   //  placing a mutation at a leaf only the sample at the leaf has the mutation
//	for(int sample=0; sample<m; sample++){
//		if(sample != leaf){
//			score += logScores[obsMutProfiles[sample][mut]][0];  // all other samples do not have it
//		}
//	}
//	return score;
//}


/* computes the best score for placing a mutation in a given binary tree */
//double binTreeMutScore(int** obsMutProfiles, int mut, int m, double ** logScores, bool** ancMatrix){
//
//	int nodeCount = (2*m)-1;
//	double bestPlacementScore = binTreeRootScore(obsMutProfiles, mut, m, logScores);
//	//print_boolMatrix(bool** array, int n, int m);
//	for(int p=0; p<nodeCount-1; p++){                           // try all possible placements (nodes in the mutation tree)
//
//		double score = 0.0;                   // score for placing mutation at a specific node
//		for(int sample=0; sample<m; sample++){
//			//cout << p << " " << sample << "\n";
//			if(ancMatrix[p][sample] == 1){
//				score += logScores[obsMutProfiles[sample][mut]][1]; // sample should have the mutation
//			}
//			else{
//				score += logScores[obsMutProfiles[sample][mut]][0]; // sample should not have the mutation
//			}
//		}
//		bestPlacementScore = max(bestPlacementScore, score);
//	}
//	return bestPlacementScore;
//}
