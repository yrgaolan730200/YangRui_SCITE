/*
 * treelist.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 */

/*	treelist.cpp

1. 简介： 
	管理和维护 MCMC 搜索过程中发现的最优树（最高似然得分树）列表及其关联参数。

2. 背景说明： 
	在 SCITE 算法的 MCMC 采样过程中，模型需要记录所有具有最高似然得分的等价树结构。该文件提供了更新、去重
	和存储这些最优树及其错误率参数（beta）的逻辑，对应于算法在搜索空间中保留全局最优解的步骤。

3.主要内容：
	updateTreeList()：根据当前得分决定是重置最优树列表还是向列表中添加新的等效树。 
	resetTreeList()：当发现得分更高的新树时，清空当前列表并插入该新树。 emptyVectorFast()：释放 vector 中存储的所有 treeBeta 结构体的动态内存并清空容器。 emptyTreeList()：释放 vector 中存储的所有整型数组（树结构）内存。 createNewTreeListElement()：深拷贝树结构并结合 beta 参数创建一个新的 treeBeta 结构体。 isDuplicateTreeFast()：检查给定的树结构是否已经存在于当前的最优树列表中。

4.输入：
	bestTrees：存储最优树及其 beta 参数的 struct treeBeta 向量引用。 
	currTreeParentVec：当前 MCMC 步生成的树父节点向量。 currScore：当前树的对数似然得分。 
	bestScore：目前为止发现的最高对数似然得分。 beta：与当前树关联的假阴性错误率参数。

5.输出：
	bestTrees：更新后的最优树集合。 
	返回值：大部分函数为 void，通过引用或指针直接修改列表。

6.注意事项：
	算法假设：
		树结构通过父节点向量（Parent Vector）表示，其中包含 n 个突变节点；判定重复树时仅比对父节点向量的数值一致性。 
	使用时的重要限制或潜在问题：
		函数内部涉及大量的动态内存分配（new/delete），需确保在列表重置或程序结束时调用相应的清理函数以防止内存泄漏；
		去重检查 isDuplicateTreeFast 采用线性搜索，当等价树数量极大时可能存在性能瓶颈。
*/

#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <sstream>
#include <math.h>
#include <queue>
#include "matrices.h"
#include "treelist.h"
#include "rand.h"

using namespace std;


void updateTreeList(vector<struct treeBeta>& bestTrees, int* currTreeParentVec, int n, double currScore, double bestScore, double beta){

	if(currScore > bestScore){
		//cout << "tree list of size " << bestTrees.size() << " emptied\n";
		resetTreeList(bestTrees, currTreeParentVec, n, beta);                              // empty the list of best trees and insert current tree

	}
	else if (currScore == bestScore){
		if(!isDuplicateTreeFast(bestTrees, currTreeParentVec, n)){               // if the same tree was not previously found
			treeBeta newElem = createNewTreeListElement(currTreeParentVec, n, beta);
			bestTrees.push_back(newElem);        // add it to list
		}
	}
}


/* removes all elements from the vector and inserts the new best tree */
void resetTreeList(vector<struct treeBeta>& bestTrees, int* newBestTree, int n, double beta){
	emptyVectorFast(bestTrees, n);                                         // empty the list of best trees
	treeBeta newElem = createNewTreeListElement(newBestTree, n, beta);
	bestTrees.push_back(newElem);                // current tree is now the only best tree
}


/* removes all elements from the vector */
void emptyVectorFast(std::vector<struct treeBeta>& optimalTrees, int n){
    for(int i=0; i<optimalTrees.size(); i++){
    	delete [] optimalTrees[i].tree;
	}
    optimalTrees.clear();
}

/* removes all elements from the vector */
void emptyTreeList(std::vector<int*>& optimalTrees, int n){
    for(int i=0; i<optimalTrees.size(); i++){
    	delete [] optimalTrees[i];
	}
    optimalTrees.clear();
}

/* creates a new tree/beta combination */
struct treeBeta createNewTreeListElement(int* tree, int n, double beta){
	treeBeta newElem;
	newElem.tree = deepCopy_intArray(tree, n);
	newElem.beta = beta;
	return newElem;
}

/* returns true if the same tree was found before */
bool isDuplicateTreeFast(std::vector<struct treeBeta> &optimalTrees, int* newTree, int n){
    for(int k=0; k<optimalTrees.size(); k++){
      bool same = true;
      for(int i=0; i<n; i++){
    	  if(newTree[i] != optimalTrees[k].tree[i]){
              same = false;
              break;
          }
      }
      if(same == true){
        return true;
      }
    }
    return false;
}







///* returns true if the tree has a branching point */
//bool hasBranching(int* parents, int n){
//	bool* isParent = new bool[n+1];
//	for(int i=0; i<n; i++){
//		isParent[i]=false;
//	}
//	for(int i=0; i<n; i++){
//		if(isParent[parents[i]]==true){
//			return true;
//		}
//		isParent[parents[i]] = true;
//	}
//	return false;
//}
//
//void foundBranchingTree(std::vector<bool**> treeList, int n){
//	bool isBranching = false;
//	vector<int> branchingTrees;
//	for(int i=0; i<treeList.size(); i++){
//		int* parVector = ancMatrixToParVector(treeList[i], n);
//		if(hasBranching(parVector, n)){
//			branchingTrees.push_back(i);
//			isBranching = true;
//		}
//	}
//	if(isBranching==true){
//		cout << branchingTrees.size() << " out of " << treeList.size() << "trees are branching: ";
//		for (int i=0; i<branchingTrees.size(); i++){
//		    cout << branchingTrees[i] << " ";
//		}
//		cout << "\n";
//	}
//	else{
//		printf("no branching trees in list\n");
//	}
//}

//bool isDuplicateTree(std::vector<bool**> &optimalTrees, bool** newTree, int n){
//    for(int k=0; k<optimalTrees.size(); k++){
//      bool same = true;
//      for(int i=0; i<n && same==true; i++){
//        for(int j=0; j<n && same==true; j++){
//            if(newTree[i][j] != optimalTrees[k][i][j]){
//              same = false;
//            }
//        }
//      }
//      if(same == true){
//        return true;
//      }
//    }
//    return false;
//}

//double updateListOfBestTrees(double currScore, double bestScore, bool**& propAncMatrix, std::vector<bool**> &optimalTrees, int n){
//
//    if(currScore > bestScore){                                    // a tree with better score is found
//        bestScore = currScore;
//        emptyVector(optimalTrees, n);     // empty outdated list of best trees
//        optimalTrees.push_back(propAncMatrix);	               // add new optimal tree to list
//  	}
//  	else if(currScore == bestScore){                            // another instance with current best score is found
//  		  optimalTrees.push_back(propAncMatrix);	                    // add the new optimal tree to the list
//  	}
//    else{
//        free_boolMatrix(propAncMatrix);
//    }
//	  return bestScore;
//    return 0;
//}

///* empty vector */
//void emptyVector(std::vector<bool**> &optimalTrees, int n){
//    for(int i=optimalTrees.size()-1; i>=0; i--){
//    	free_boolMatrix(optimalTrees.at(i));
//    	optimalTrees.at(i) = NULL;
//	  }
//	  optimalTrees.clear();
//}
