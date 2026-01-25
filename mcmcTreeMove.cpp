/*
 * mcmcTreeMove.cpp
 *
 *  Created on: Mar 15, 2016
 *      Author: jahnka
 */

 /*mcmcTreeMove.cpp

 1. 简介：
    实现突变树在 MCMC 采样过程中的各种拓扑结构变动（Moves）策略，包括剪枝重接、节点交换及子树交换。

 2.背景说明：
    该文件是 SCITE 等生物信息学推断工具的核心组件。它定义了如何在突变进化树（Mutation Tree）的空间内进行随机搜索。
    这些变动算子保证了 MCMC 链在树空间的遍历性，是计算最大后验概率树或进行贝叶斯推断的关键步骤。

 3.主要内容：
    proposeNewTree()：根据给定的概率分布随机选择并执行一种变动算子（剪枝、交换标签或交换子树）来生成新候选树
    choseParent()：从候选父节点集中（包括根节点）随机选取一个节点作为目标挂载点
    getNewParentVecFast()：通过修改特定节点的父节点快速生成新的父节点向量
    getNewParentVec_SwapFast()：通过交换两个节点的标签及相应的父子关系快速生成新向量
    reorderToStartWithDescendant()：调整两个节点的顺序，确保在同一谱系时先处理后代节点，以维持树的结构合法性
    getNewAncMatrix_Swap()：通过行列交换更新节点的祖先矩阵
    getNewAncMatrix()：在剪枝重接操作后，根据新的拓扑关系增量更新祖先矩阵
    
 4.输入：
    moveProbs：各种变动类型（Prune & Re-attach, Swap Labels, Swap Subtrees）的选择概率向量
    n：树中节点的数量（通常代表突变位点数）
    currTreeAncMatrix：当前树的祖先矩阵（$n \times n$）
    currTreeParentVec：当前树的父节点向量表示nbhcorrection：引用参数，用于存储并返回变动导致的邻域修正因子（Neighborhood Correction）

 5.输出：
    propTreeParVec：生成并返回新候选树的父节点向量指针（若变动失败则返回 NULL）
    
 6.注意事项：
    算法假设：
        假设生成的结构始终为有根树（Rooted Tree），根节点通常表示为 $n$ 或 $n+1$。

    使用时的重要限制或潜在问题：
        在交换位于同一谱系（Lineage）的子树时，必须引入邻域修正因子以修正转移概率的不对称性，
        否则会破坏 MCMC 的细致平衡条件（Detailed Balance）。代码中通过 nbhcorrection 实现此项修正。
*/
 

#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include "matrices.h"
#include "treelist.h"
#include "trees.h"
#include "scoreTree.h"
#include "rand.h"
#include "limits.h"
#include "output.h"
#include "mcmcTreeMove.h"

using namespace std;

int* proposeNewTree(vector<double> moveProbs, int n, bool** currTreeAncMatrix, int* currTreeParentVec, double& nbhcorrection){

	int* propTreeParVec = NULL;
	int movetype = sampleRandomMove(moveProbs);      // pick the move type
	//cout << "move type: " << movetype << "\n";
	nbhcorrection = 1;                               // reset the neighbourhood correction

	if(movetype==1){       /* prune and re-attach */
		//cout << "move type is prune and reattach\n";
		int nodeToMove = pickRandomNumber(n);   // pick a node to move with its subtree
		std::vector<int> possibleparents = getNonDescendants(currTreeAncMatrix, nodeToMove, n);      // possible attachment points
		int newParent = choseParent(possibleparents, n);                                             // randomly pick a new parent among available nodes, root (n+1) is also possible parent
		propTreeParVec = getNewParentVecFast(currTreeParentVec, nodeToMove, newParent, n);           // create new parent vector
	}
    else if(movetype==2){   /* swap two node labels  */
    	//cout << "move type is swap node labels\n";
    	int* nodestoswap = sampleTwoElementsWithoutReplacement(n);
        propTreeParVec    = getNewParentVec_SwapFast(currTreeParentVec, nodestoswap[0], nodestoswap[1], n);
        delete [] nodestoswap;
    }
    else if(movetype==3){    /*  swap two subtrees  */
	   //cout << "move type is swap subtrees\n";
        int* nodestoswap = sampleTwoElementsWithoutReplacement(n);                    // pick the node that will be swapped
        nodestoswap = reorderToStartWithDescendant(nodestoswap, currTreeAncMatrix);   // make sure we move the descendant first (in case nodes are in same lineage)
        int nodeToMove = nodestoswap[0];                                              // now we move the first node chosen and its descendants
        int nextnodeToMove = nodestoswap[1];                                          // next we need to move the second node chosen and its descendants
        delete [] nodestoswap;

        if(currTreeAncMatrix[nextnodeToMove][nodeToMove]==0){                      // the nodes are in different lineages -- simple case

        	propTreeParVec = deepCopy_intArray(currTreeParentVec, n);         // deep copy of parent matrix to keep old one
        	propTreeParVec[nodeToMove] = currTreeParentVec[nextnodeToMove];        // and exchange the parents of the nodes
        	propTreeParVec[nextnodeToMove] = currTreeParentVec[nodeToMove];
        }
        else{                                                                      // the nodes are in the same lineage -- need to avoid cycles in the tree
        	propTreeParVec = deepCopy_intArray(currTreeParentVec, n);         // deep copy of parent vector to keep old one
        	propTreeParVec[nodeToMove] = currTreeParentVec[nextnodeToMove];        // lower node is attached to the parent of the upper node
        	std::vector<int> descendants     = getDescendants(currTreeAncMatrix, nodeToMove, n);   // all nodes in the subtree of the lower node
        	bool** propTreeAncMatrix = parentVector2ancMatrix(propTreeParVec, n);
        	std::vector<int> nextdescendants = getDescendants(propTreeAncMatrix, nextnodeToMove, n);
        	free_boolMatrix(propTreeAncMatrix);
        	propTreeParVec[nextnodeToMove] = descendants[pickRandomNumber(descendants.size())];  // node to move is attached to a node chosen uniformly from the descendants of the first node
        	nbhcorrection = 1.0*descendants.size()/nextdescendants.size(); // neighbourhood correction needed for MCMC convergence, but not important for simulated annealing
        }
    }
    return propTreeParVec;
}


/* picks a parent randomly from the set of possible parents, this set includes the root (n+1) */
int choseParent(std::vector<int> &possibleParents, int root){
	possibleParents.push_back(root);                           // add root, as it is also possible attachement point
    int chosenParentPos = pickRandomNumber(possibleParents.size());  // choose where to append the subtree
    int newParent = possibleParents[chosenParentPos];
	possibleParents.pop_back();    // remove root from list of possible parents as it is treated as special case later on
	return newParent;
}



/* creates the new parent vector after pruning and reattaching subtree */
int* getNewParentVecFast(int* currTreeParentVec, int nodeToMove, int newParent, int n){
	int* propTreeParVec = deepCopy_intArray(currTreeParentVec, n);        // deep copy of parent matrix to keep old one
    propTreeParVec[nodeToMove] = newParent;                       // only the parent of the moved node changes
    return propTreeParVec;
}


/* creates the new parent vector after swapping two nodes */
int* getNewParentVec_SwapFast(int* currTreeParentVec, int first, int second, int n){

	int* propTreeParVec = deepCopy_intArray(currTreeParentVec, n);
    for(int i=0; i<n; i++){           // set vector of proposed parents
      if(propTreeParVec[i] == first && i!=second){
      	propTreeParVec[i] = second;              // update entries with swapped parents
      }
      else if(propTreeParVec[i] == second && i!=first){
      	propTreeParVec[i] = first;
      }
    }
    int temp = propTreeParVec[first];
    propTreeParVec[first] = propTreeParVec[second];  // update parents of swapped nodes
    propTreeParVec[second] = temp;
    if(propTreeParVec[first]==first){propTreeParVec[first]=second;}    // this is needed to ensure that tree is connected, the above fails
    if(propTreeParVec[second]==second){propTreeParVec[second]=first;}  // if first is parent of second, or vice versa
    return propTreeParVec;
}

/* re-orders the nodes so that the descendant is first in case the nodes are in the same lineage */
int* reorderToStartWithDescendant(int* nodestoswap, bool** currTreeAncMatrix){
    if(currTreeAncMatrix[nodestoswap[0]][nodestoswap[1]]==true){  // make sure we move the descendent first
        int temp = nodestoswap[0];
        nodestoswap[0] = nodestoswap[1];
        nodestoswap[1] = temp;
    }
    return nodestoswap;
}


/* creates the new parent vector after swapping two nodes */
int* getNewParentVec_Swap(int* currTreeParentVec, int first, int second, int n, int* propTreeParVec){
    for(int i=0; i<n; i++){           // set vector of proposed parents
      if(propTreeParVec[i] == first && i!=second){
      	propTreeParVec[i] = second;              // update entries with swapped parents
      }
      else if(propTreeParVec[i] == second && i!=first){
      	propTreeParVec[i] = first;
      }
    }

    int temp = propTreeParVec[first];
    propTreeParVec[first] = propTreeParVec[second];  // update parents of swapped nodes
    propTreeParVec[second] = temp;
    if(propTreeParVec[first]==first){propTreeParVec[first]=second;}    // this is needed to ensure that tree is connected, the above fails
    if(propTreeParVec[second]==second){propTreeParVec[second]=first;}  // if first is parent of second, or vice versa
    return propTreeParVec;
}


/* creates the new ancestor matrix after swapping two nodes, old matrix is kept */
bool** getNewAncMatrix_Swap(bool** currTreeAncMatrix, int first, int second, int n, bool** propTreeAncMatrix){

    for(int i=0; i<n; i++){                                       // swap columns
    	bool temp = propTreeAncMatrix[i][first];
      	propTreeAncMatrix[i][first] = propTreeAncMatrix[i][second];
      	propTreeAncMatrix[i][second] = temp;
      }
      for(int i=0; i<n; i++){                                // swap rows
      	bool temp = propTreeAncMatrix[first][i];
      	propTreeAncMatrix[first][i] = propTreeAncMatrix[second][i];
      	propTreeAncMatrix[second][i] = temp;
      }
     return propTreeAncMatrix;
 }


/* creates the new parent vector after pruning and reattaching subtree */
int* getNewParentVec(int* currTreeParentVec, int nodeToMove, int newParent, int n, int *propTreeParVec){
    propTreeParVec[nodeToMove] = newParent;                       // only the parent of the moved node changes
    return propTreeParVec;
}


/* creates the new ancestor matrix after pruning and reattaching subtree, old matrix is kept */
bool** getNewAncMatrix(bool** currTreeAncMatrix, int newParent, std::vector<int> descendants, std::vector<int> possibleParents, int n, bool** propTreeAncMatrix){

    if(newParent<n){    // replace the non-descendants of the node and its descendants by the ancestors of the new parent
		for(int i=0; i<possibleParents.size(); i++){
  			for(int j=0; j<descendants.size(); j++){
  				propTreeAncMatrix[possibleParents[i]][descendants[j]] = currTreeAncMatrix[possibleParents[i]][newParent];
  			}
  		}
    }
    else
    {     // if we attach to the root, then they have no further ancestors
        for(int i=0; i<possibleParents.size(); i++){
  			for(int j=0; j<descendants.size(); j++){
        		propTreeAncMatrix[possibleParents[i]][descendants[j]] = 0;
        	}
        }
    }
    return propTreeAncMatrix;
}


//int proposeNextTree(std::vector<double> moveProbs, int n, bool** currTreeAncMatrix, int* currTreeParentVec, int*& propTreeParVec, bool**& propTreeAncMatrix){
//
//	int movetype = sampleRandomMove(moveProbs);   // pick the move type
//	int nbhcorrection = 1;               // reset the neighbourhood correction
//	if(movetype==1){       /* prune and re-attach */
//		//cout << "move type is prune and reattach\n";
//		int nodeToMove = pickRandomNumber(n);   // pick a node to move with its subtree
//
//		//  descendants are those in row nodeToMove with a 1, possible parents are the others
//		std::vector<int> descendants  = getDescendants(currTreeAncMatrix, nodeToMove, n);                   // includes the node itself
//		std::vector<int> possibleparents = getNonDescendants(currTreeAncMatrix, nodeToMove, n);             // possible attachment points
//
//		int newParent = choseParent(possibleparents, n);                                    // randomly pick a new parent among available nodes, root (n+1) is also possible parent
//		propTreeAncMatrix = getNewAncMatrix(currTreeAncMatrix, newParent, descendants, possibleparents, n, propTreeAncMatrix);   // create new ancestor matrix, TODO: avoid deepCopy of old matrix
//		propTreeParVec    = getNewParentVec(currTreeParentVec, nodeToMove, newParent, n, propTreeParVec);                     // create new parent vector
//	}
//    else if(movetype==2){   /* swap two node labels  */
//    //	cout << "move type is swap node labels\n";
//    	int* nodestoswap = sampleTwoElementsWithoutReplacement(n);
//        propTreeAncMatrix = getNewAncMatrix_Swap(currTreeAncMatrix, nodestoswap[0], nodestoswap[1], n, propTreeAncMatrix);  // create ancestor matrix after swapping
//        propTreeParVec    = getNewParentVec_Swap(currTreeParentVec, nodestoswap[0], nodestoswap[1], n, propTreeParVec);
//        delete [] nodestoswap;
//    }
//   else{       // if(movetype==3){    /*  swap two subtrees  */
//	  // cout << "move type is swap subtrees\n";
//        int* nodestoswap = sampleTwoElementsWithoutReplacement(n);                   // pick the nodes whose incoming edges will be cut and swapped around
//        nodestoswap = reorderToStartWithDescendant(nodestoswap, currTreeAncMatrix);      // make sure we move the descendant first (in case nodes are in same lineage)
//        int nodeToMove = nodestoswap[0];                                              // now we move the first node chosen and its descendants
//        //cout << "picked nodes to prune tree: " << nodestoswap[0]+1 << " and " << nodestoswap[1]+1 << "\n";
//        std::vector<int> descendants  = getDescendants(currTreeAncMatrix, nodeToMove, n);   // includes the node itself
//        std::vector<int> possibleparents = getNonDescendants(currTreeAncMatrix, nodeToMove, n);  // possible attachment points
//        int newParent = currTreeParentVec[nodestoswap[1]];                                           // the new parent is the parent of the other node
//        //cout << "new parent: " << newParent+1 << " of " <<  nodestoswap[0]+1 << "\n";
//        bool** propTreeAncMatrixTemp = getNewAncMatrix(currTreeAncMatrix, newParent, descendants, possibleparents, n, propTreeAncMatrix);
//        propTreeParVec = getNewParentVec(currTreeParentVec, nodeToMove, newParent, n, propTreeParVec);
//       // print_intArray(propTreeParVec, n);
//        int nextnodeToMove = nodestoswap[1];  // next we need to move the second node chosen and its descendants
//        delete [] nodestoswap;
//        std::vector<int> nextdescendants     = getDescendants(propTreeAncMatrixTemp, nextnodeToMove, n);  // includes the node itself
//        std::vector<int> nextpossibleparents = getNonDescendants(propTreeAncMatrixTemp, nextnodeToMove, n);
//
//        if(currTreeAncMatrix[nextnodeToMove][nodeToMove]==1){             // if the two nodes used to be in the same lineage (first descendant of second) , the second
//        	newParent = descendants[pickRandomNumber(descendants.size())];  // node to move is attached to a node chosen uniformly from the descendants of the first node
//        	nbhcorrection = descendants.size()/nextdescendants.size(); // neighbourhood correction needed for MCMC convergence, but not important for simulated annealing
//        }
//        else  // if the nodes used to be in different lineages, just attach to the second node to initial parent of the first node
//        {
//        	newParent = currTreeParentVec[nodeToMove];
//        }
//
//        propTreeAncMatrix = getNewAncMatrix(propTreeAncMatrixTemp, newParent, nextdescendants, nextpossibleparents, n, deepCopy_boolMatrix(propTreeAncMatrixTemp, n, n));
//        free_boolMatrix(propTreeAncMatrixTemp);
//        propTreeParVec[nextnodeToMove] = newParent; // update the parent vector of the moved node
//   }
//	return nbhcorrection;
//}
