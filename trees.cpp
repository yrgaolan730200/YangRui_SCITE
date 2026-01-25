/*
 * trees.cpp
 *
 *  Created on: Oct 12, 2015
 *      Author: jahnka
 */
/*	trees.cpp
1. 简介：
	实现树结构的各种表示形式（父节点向量、祖先矩阵、Prüfer 序列、子节点列表）之间的转换及树的拓扑操作。

2.背景说明：
	该文件是 SCITE 算法中处理进化树拓扑的核心工具库。它负责在 MCMC 采样过程中生成初始随机树、维护树的层次关系，
	并为似然得分计算和最终的 Newick 格式导出提供底层算法支持，主要对应论文中关于搜索空间拓扑表示和变换的部分。

3. 主要内容：
	getDescendants()：通过祖先矩阵获取指定节点的所有后代节点集合。
	getNonDescendants()：获取所有不是指定节点后代的节点。
	countBranches()：通过统计子节点列表中的叶子节点数量来计算树的支路数。
	getChildListFromParentVector()：将父节点向量转换为邻接表格式的子节点列表。
	getNewickCode()：递归生成符合 Newick 标准的树描述字符串，用于系统发育分析。
	getBreadthFirstTraversal()：计算树的广度优先遍历序列（BFT），常用于优化评分时的顺序访问。
	parentVector2ancMatrix()：将父节点向量转换为完整的祖先关系矩阵（包含自环）。
	prueferCode2parentVector()：使用线性时间的 Prüfer 序列算法重构对应的唯一树父节点向量。
	getRandParentVec()：生成一个包含 n 个突变位点且以第 n 号节点为根的随机树父节点向量。
	starTreeVec() / starTreeMatrix()：生成星形树（Star Tree）的向量或矩阵表示，通常用于搜索的初始状态。

4.	输入：
	parent：存储节点父代索引的整型数组。
	ancMatrix：表示节点间祖先关系的布尔型矩阵。
	code/codeLength：Prüfer 编码序列及其长度。
	n：突变位点或节点的数量。
	root：作为递归起点或树根的节点索引。

5.	输出：
	childList：嵌套的 vector 结构，表示每个节点的子节点集合。
	ancMatrix：转换后的 $n \times n$ 祖先矩阵。
	parentVector：重构或生成的树父节点向量指针。
	Newick string：树的 Newick 格式字符串。

6. 注意事项：
	算法假设：
		假定树是连通且无环的；在 prueferCode2parentVector 中，假定根节点索引始终为 codeLength + 1。
	使用时的重要限制或潜在问题：
		返回的数组（如 bft, parent）通常是在堆上分配的，调用者需负责内存释放；getNewickCode 必须从真实的根节点开始递归，否则生成的字符串不完整。
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
#include "trees.h"
#include "rand.h"
#include "output.h"

using namespace std;


/* returns all nodes that are descendants of the given node */
/* note: ancMatrix is 1 at [i,j] if i is an ancestor of j in the tree */
std::vector<int> getDescendants(bool** ancMatrix, int node, int n){
  std::vector<int> descendants;
  for(int i=0; i<n; i++){
  	if(ancMatrix[node][i]==true){
			descendants.push_back(i);
		}
	}
	return descendants;
}

/* returns all nodes that are not descendants of the given node */
/* i.e. ancestors and nodes in a different branch of the tree   */
/* note: ancMatrix is 0 at [i,j] if i is not an ancestor of j in the tree */
std::vector<int> getNonDescendants(bool**& ancMatrix, int node, int n){
	std::vector<int> ancestors;
	for(int i=0; i<n; i++){
		if(ancMatrix[node][i]==false){
			ancestors.push_back(i);
		}
	}
	return ancestors;
}

/* counts the number of branches in a tree, this is the same as the number of leafs in the tree */
int countBranches(int* parents, int length){
	int count = 0;
	vector<vector<int> > childList = getChildListFromParentVector(parents, length);
	for(int i=0; i<childList.size(); i++){
		if(childList.at(i).size()==0){ count++; }
	}
	for(int i=0; i<childList.size(); i++){
		childList[i].clear();
	}
	childList.clear();
	return count;
}

/* converts a parent vector to the list of children */
vector<vector<int> > getChildListFromParentVector(int* parents, int n){

	vector<vector<int> > childList(n+1);
	for(int i=0; i<n; i++){
		childList.at(parents[i]).push_back(i);
	}
	return childList;
}

void deleteChildLists(vector<vector<int> > &childLists){
	for(int i=0; i<childLists.size(); i++){
		childLists[i].clear();
	}
	childLists.clear();
}

/* converts a tree given as lists of children to the Newick tree format */
/* Note: This works only if the recursion is started with the root node which is n+1 */
string getNewickCode(vector<vector<int> > list, int root){
	stringstream newick;
	vector<int> rootChilds = list.at(root);
	if(!rootChilds.empty()){
		newick << "(";
		bool first = true;
		for(int i=0; i<rootChilds.size(); i++){
			if(!first){
				newick << ",";
			}
			first = false;
			newick << getNewickCode(list, rootChilds.at(i));
		}
		newick << ")";
	}
	newick << root+1;
	return newick.str();
}



/*  computes a breadth first traversal of a tree from the parent vector  */
int* getBreadthFirstTraversal(int* parent, int n){

	vector<vector<int> > childLists = getChildListFromParentVector(parent, n);
	int* bft = new int[n+1];
	bft[0] = n;
	int k = 1;

	for(int i=0; i<n+1; i++){
		for(int j=0; j<childLists[bft[i]].size(); j++){
			bft[k++] = childLists[bft[i]][j];
		}
	}
	for(int i=0; i<childLists.size(); i++){
		childLists[i].clear();
	}
	childLists.clear();
	return bft;
}

int* reverse(int* array, int length){
	int temp;

	for (int i = 0; i < length/2; ++i) {
		temp = array[length-i-1];
		array[length-i-1] = array[i];
		array[i] = temp;
	}
	return array;
}



/* transforms a parent vector to an ancestor matrix*/
bool** parentVector2ancMatrix(int* parent, int n){
	bool** ancMatrix = init_boolMatrix(n, n, false);
	int root = n;
	for(int i=0; i<n; i++){
		int anc = i;
		int its =0;
		while(anc < root){                              // if the ancestor is the root node, it is not represented in the adjacency matrix
			if(parent[anc]<n){
				ancMatrix[parent[anc]][i] = true;
			}

			anc = parent[anc];
			its++;
		}
	}
	for(int i=0; i<n; i++){
		ancMatrix[i][i] = true;
	}
	return ancMatrix;
}

/* given a Pruefer code, compute the corresponding parent vector */
int* prueferCode2parentVector(int* code, int codeLength){
	int nodeCount = codeLength+1;
	int* parent = new int[nodeCount];
	//print_intArray(code, codeLength);
	int* lastOcc = getLastOcc(code, codeLength);    // node id -> index of last occ in code, -1 if no occurrence or if id=root
	bool* queue = getInitialQueue(code, codeLength);  // queue[node]=true if all children have been attached to this node, or if it is leaf
	int queueCutter = -1;    // this is used for a node that has been passed by the "queue" before all children have been attached
	int next = getNextInQueue(queue, 0, codeLength+1);

	for(int i=0; i<codeLength; i++){               // add new edge to tree from smallest node with all children attached to its parent
		if(queueCutter >=0){
			parent[queueCutter] = code[i];         // this node is queueCutter if the queue has already passed this node
			//cout << queueCutter << " -> " << code[i] << "\n";
			queueCutter = -1;
		}
		else{
			parent[next] = code[i];                               // use the next smallest node in the queue, otherwise
			//cout << next << " -> " << code[i] << "\n";
			next = getNextInQueue(queue, next+1, codeLength+1);     // find next smallest element in the queue
		}

		if(lastOcc[code[i]]==i){                               // an element is added to the queue, or we have a new queueCutter
			updateQueue(code[i], queue, next);
			queueCutter = updateQueueCutter(code[i], queue, next);
		}
	}
	if(queueCutter>=0){
		parent[queueCutter] = nodeCount;
		//cout << queueCutter << " -> " << nodeCount << "\n";
	}
	else{
		parent[next] = nodeCount;
		//cout << next << " -> " << nodeCount << "\n";
	}

	delete [] lastOcc;
	delete [] queue;
	//print_intArray(parent, codeLength+1);
	//getGraphVizFileContentNumbers(parent, codeLength+1);
	return parent;
}

bool* getInitialQueue(int* code, int codeLength){
	//cout << "code Length: " << codeLength << "\n";
	int queueLength = codeLength+2;
	//cout << "queueLength: " << queueLength << "\n";
	bool* queue = init_boolArray(queueLength, true);

	for(int i=0; i<codeLength; i++){
		queue[code[i]] = false;
	}
	return queue;
}


void updateQueue(int node, bool* queue, int next){

	if(node>=next){                //  add new node to queue
		queue[node] = true;
	}
}

int updateQueueCutter(int node, bool* queue, int next){
	if(node>=next){
		return -1;         // new node can be added to the queue
	}
	else{
		return node;         // new node needs to cut the queue, as it has already passed it
	}
}


int* getLastOcc(int* code, int codeLength){
	int* lastOcc = init_intArray(codeLength+2, -1);
	int root = codeLength+1;
	for(int i=0; i<codeLength; i++){
		if(code[i] != root){
			lastOcc[code[i]] = i;
		}
	}
	return lastOcc;
}

int getNextInQueue(bool* queue, int pos, int length){
	for(int i=pos; i<length; i++){
		if(queue[i]==true){
			return i;
		}
	}
	//cout << "No node left in queue. Possibly a cycle?";
	return length;
}

/* creates a random parent vector for nodes 0, .., n with node n as root*/
int* getRandParentVec(int n){
	int* randCode = getRandTreeCode(n);
	int* randParent = prueferCode2parentVector(randCode, n-1);
	delete [] randCode;
	return randParent;
}



/* creates the parent vector for a star tree with node n as center and 0,...,n-1 as leafs */
int* starTreeVec(int n){
	int* starTreeVec = new int[n];
	for(int i=0;i<n;i++){
		starTreeVec[i] = n;
	}
	return starTreeVec;
}

/* creates the ancestor matrix for the same tree */
bool** starTreeMatrix(int n){

  bool** starTreeMatrix = init_boolMatrix(n, n, false);
  for(int i=0;i<n;i++){
		starTreeMatrix[i][i] = true;
	}
	return starTreeMatrix;
}
