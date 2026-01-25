/*
 * mcmcBinTreeMove.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: jahnka
 */

 /*	mcmcBinTreeMove.cpp

1. 简介： 
	实现针对二叉树结构的 MCMC 提案（Moves）逻辑，通过剪枝重接或交换叶节点标签来生成新的候选树。

2. 背景说明： 
	该文件主要用于 SCITE 算法在“转置”模式（-transpose）下的运行。在此模式下，算法推断的是样本（细胞）间的二叉演化树，而非突变树。
	该代码对应于肿瘤系统发育树重构中在二叉树拓扑空间进行随机搜索的步骤。

3. 主要内容： 
	proposeNextBinTree()：根据设定的概率随机选择并执行一种二叉树变动操作，生成新树的父节点向量。 
	pickNodeToMove()：在树中随机选择一个可移动的合法节点（要求其父节点非根节点）。 
	getSibling()：获取二叉树中指定节点的唯一兄弟节点，用于在剪枝操作后维护树的二叉性质。

4. 输入： 
	moveProbs：各种变动类型（剪枝重接、标签交换）的概率分布向量。 
	m：叶节点（细胞）的数量。 
	currTreeParVec：当前二叉树的父节点向量表示。 
	currTreeAncMatrix：当前二叉树的祖先矩阵，用于快速判断节点间的拓扑关系。

5. 输出： 
	propTreeParVec：执行变动操作后产生的新树的父节点向量指针。

6. 注意事项： 
	算法假设：
		假设所处理的树是严格的二叉树（除了根节点外，每个非叶节点都有且仅有两个子节点）。 
	使用时的重要限制或潜在问题：
		剪枝重接（SPR）操作会将节点与其父节点作为一个整体移动，以保持二叉结构；标签交换目前仅限于前 m 个节点（即叶节点/细胞）。
*/
#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <iostream>
#include "matrices.h"
#include "trees.h"
#include "rand.h"
#include "mcmcBinTreeMove.h"
#include "output.h"
using namespace std;


/* proposes a new binary tree by a single move from the current binary tree based on the move probabilities */
/* the old tree is kept as currTree, the new one is stored as propTreeParVec */
int* proposeNextBinTree(std::vector<double> moveProbs, int m, int* currTreeParVec, bool** currTreeAncMatrix){

	int movetype = sampleRandomMove(moveProbs);      // pick the move type according to move probabilities
	int parVecLength = (2*m)-2;               // 2m-1 nodes, but the root has no parent

	//cout << "move prob 0: " << moveProbs[0] << "\n";
	//cout << "move prob 1: " << moveProbs[1] << "\n";
	//cout << "move prob 2: " << moveProbs[2] << "\n";
	vector<vector<int> >childLists = getChildListFromParentVector(currTreeParVec, parVecLength);
	int* propTreeParVec  = deepCopy_intArray(currTreeParVec, parVecLength);

	if(movetype==1){                                                       /* type 1: prune and re-attach */
		//cout << "move type is prune and re-attach in binary tree\n";
		int v = pickNodeToMove(currTreeParVec, parVecLength);
		int p = currTreeParVec[v];
		int sib = getSibling(v, currTreeParVec, childLists);             // get the sibling of node v and attach it to the
		propTreeParVec[sib] = currTreeParVec[p];                         // grandparent of v, as the parent of v is moved along with v

		std::vector<int> possibleSiblings = getNonDescendants(currTreeAncMatrix, p, parVecLength);    // get list of possible new siblings of v

		if(possibleSiblings.size()==0){
			cerr << "Error: No new sibling found for node " << v << " for move type 1 in binary tree.\n"; // Should never occur. Something wrong with the tree.
			printGraphVizFile(currTreeParVec, parVecLength);
		}

		int newSibling = possibleSiblings[pickRandomNumber(possibleSiblings.size())]; // pick a new sibling from remaining tree (root can not be a sibling)
		propTreeParVec[newSibling] = p;                                               // make the new sibling a child of v's parent
		propTreeParVec[p] = currTreeParVec[newSibling];                            // make the parent of v the child of the new sibling's former parent
	}
    else{                                                                 /* type 2: swap two node labels  */
    	//cout << "move type is swap node labels in binary tree\n";
    	int v =  rand() % m;                                            // get random leaf to swap (only the first m nodes are leafs)
    	int w =  rand() % m;                                            // get second random leaf to swap
    	propTreeParVec[v] = currTreeParVec[w];                         // and just swap parents
    	propTreeParVec[w] = currTreeParVec[v];
    }
    return propTreeParVec;
}


/* returns a node where the prune and re-attach step starts */
int pickNodeToMove(int* currTreeParentVec, int parentVectorLength){
	bool validNode = false;
	int rootId = parentVectorLength;
	while(!validNode){
		int v = pickRandomNumber(parentVectorLength);   // pick a node for the prune and re-attach step;
		if(currTreeParentVec[v]!=rootId){               // it has to be a node whose parent is not the root, as node and parent are moved together
			return v;
		}                                      // for a binary tree with more than two leafs this can not be an infinite loop
	}
}


/* returns the (unique) sibling of node v. The sibling has to exist because tree is binary */
int getSibling(int v, int* currTreeParVec, vector<vector<int> > &childLists){

	if(childLists.at(currTreeParVec[v]).at(0) != v){
		return childLists.at(currTreeParVec[v]).at(0);
	}
	else{
		return childLists.at(currTreeParVec[v]).at(1);
	}
}

