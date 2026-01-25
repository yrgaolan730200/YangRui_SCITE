/*
 * output.cpp
 *
 *  Created on: Oct 12, 2015
 *      Author: jahnka
 */

 /*	output.cpp

1. 简介： 
	处理 MCMC 搜索结果的后期处理与可视化输出，包括寻找突变的最优放置位置及生成 GraphViz 格式的树结构文件。

2. 背景说明： 
	该文件位于 SCITE 算法管线的末端，负责将抽象的 MCMC 采样结果（如父节点向量、祖先矩阵）转换为可解释的生物学结论。
	它实现了将单细胞样本挂载到进化树上的逻辑，并支持将结果导出为 .gv 文件以便进行图形化展示。

3.主要内容：
	binTreeRootScore()：计算特定突变放置在树根节点时的对数得分贡献 
	getHighestOptPlacement()：在给定样本树拓扑下，计算并返回单个突变的最优位置索引 
	getHighestOptPlacementVector()：遍历所有突变，计算每个突变在树中的最佳放置位置向量 
	getBinTreeNodeLabels()：为二叉树节点分配基因名称标签，合并同一位置的多个突变 
	getLcaWithLabel()：向上寻找最近的具有非空标签的祖先节点，用于简化图形展示 
	getGraphVizBinTree()：生成二叉样本树（Sample Tree）的 GraphViz 格式内容 
	getMutTreeGraphViz()：生成突变树（Mutation Tree）的 GraphViz 格式内容 
	writeToFile()：将生成的字符串内容写入磁盘文件 
	getGraphVizFileContentNumbers()：使用数字索引生成基础的 GraphViz 树结构 
	getGraphVizFileContentNames()：使用基因名称生成带样本挂载信息的 GraphViz 树结构 
	getBestAttachmentString()：计算样本到基因树的最佳附件关系并转换为字符串 
	attachmentPoints()：重新计算每个样本在当前树中的最佳挂载点矩阵 
	printParentVectors()：将最优树列表以父节点向量和 GraphViz 格式打印至控制台 
	printGraphVizFile()：将单一树结构的 GraphViz 文件内容打印至控制台 
	printSampleTrees()：将采样得到的树列表批量写入指定的文本文件 
	printScoreKimSimonTree()：评估并打印 Kim & Simon 方法预测树的对数似然得分

4.输入：
	obsMutProfiles/dataMatrix：观测到的单细胞突变谱矩阵 
	logScores：预计算的错误对数似然得分矩阵（依赖于假阳性/假阴性率）
	ancMatrix：表示树拓扑结构的祖先关系矩阵
	parents/currTreeParentVec：节点的父节点向量表示 
	geneNames：基因/突变位点的名称列表

5.输出：
	GraphViz string：符合 DOT 语言标准的树描述字符串 
	bestPlacements：包含各突变最优位置索引的整型数组 
	attachmentMatrix：指示样本与基因节点挂载关系的布尔矩阵

6.注意事项：
	算法假设：
		挂载逻辑假设样本应该被放置在最大化其突变观测似然的节点上；在 getHighestOptPlacement 中，若存在多个得分相同的最高位置，
		默认选择更靠近根部的节点（Highest Placement）。
	
	使用时的重要限制或潜在问题：
		GraphViz 输出主要用于中小型规模的树（如突变数 < 100），在超大规模数据下生成的图片可能难以阅读；
		某些辅助函数（如 printScoreKimSimonTree）包含硬编码的固定树结构，仅供特定的基准测试使用。
*/

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <float.h>
#include "output.h"
#include "scoreTree.h"
#include "matrices.h"

using namespace std;

/* Score contribution by a specific mutation when placed at the root, that means all samples should have it */
/* This is the same for all trees and can be precomputed */
double binTreeRootScore(int** obsMutProfiles, int mut, int m, double ** logScores){
	double score = 0.0;
	for(int sample=0; sample<m; sample++){
		score += logScores[obsMutProfiles[sample][mut]][1];
	}
	return score;
}

/* computes the best placement of a mutation, the highest one if multiple co-opt. placements exist*/
int getHighestOptPlacement(int** obsMutProfiles, int mut, int m, double ** logScores, bool** ancMatrix){

	int nodeCount = (2*m)-1;
	int bestPlacement = (2*m)-2;   // root
	double bestPlacementScore = binTreeRootScore(obsMutProfiles, mut, m, logScores);
	//cout << bestPlacementScore << " (root)\n";
	//print_boolMatrix(bool** array, int n, int m);
	for(int p=0; p<nodeCount-1; p++){                           // try all possible placements (nodes in the mutation tree)

		double score = 0.0;                   // score for placing mutation at a specific node
		for(int sample=0; sample<m; sample++){
			//cout << p << " " << sample << "\n";
			if(ancMatrix[p][sample] == 1){
				score += logScores[obsMutProfiles[sample][mut]][1]; // sample should have the mutation
			}
			else{
				score += logScores[obsMutProfiles[sample][mut]][0]; // sample should not have the mutation
			}
		}
		if(score > bestPlacementScore){
			bestPlacement = p;
			bestPlacementScore = score;
			//cout << bestPlacementScore << " (non-root)\n";
		}
		else if (score == bestPlacementScore && ancMatrix[p][bestPlacement] == true){
			bestPlacement = p;
		}
	}

	//if(bestPlacement == (2*m)-2){
	//	cout<< "best placed at root\n";
	//	getchar();
	//}
	return bestPlacement;
}

/* computes the best placement of a mutation, the highest one if multiple co-opt. placements exist*/
int* getHighestOptPlacementVector(int** obsMutProfiles, int n, int m, double ** logScores, bool** ancMatrix){
	int* bestPlacements = init_intArray(n, -1);
	for(int mut=0; mut<n; mut++){                                                               // for all mutations get
		bestPlacements[mut] = getHighestOptPlacement(obsMutProfiles, mut, m, logScores, ancMatrix);         // bestPlacementScore
	 }
	//print_intArray(bestPlacements, n);
	return bestPlacements;
}

vector<string> getBinTreeNodeLabels(int nodeCount, int* optPlacements, int n, vector<string> geneNames){
	vector<string> v;
	int count = 0;
	for(int i = 0; i < nodeCount; i++){
		v.push_back("");
	}

	for(int mut=0; mut<n; mut++){
		string toAppend;
		if(v.at(optPlacements[mut]) == ""){
			toAppend = geneNames.at(mut);
			count++;
		}
		else{
			toAppend = ", " + geneNames.at(mut);
			count++;
		}
		//cout << "        " << j << "\n";
		//cout << "                     "<< optPlacements[j] << "\n";
		v.at(optPlacements[mut]) += toAppend;
	}
	if(v.at(nodeCount-1) == ""){
		v.at(nodeCount-1) = "root";
	}
	for(int i = 0; i < nodeCount; i++){
		if(v.at(i).find(" ") != string::npos){
			v.at(i) = "\"" + v.at(i) + "\"";
		}
	}
	//cout << "added mutations " << count << "\n";
	return v;
}

/* returns the lca of a node that has a non-empty label, the root is assumed to always have a label */
int getLcaWithLabel(int node, int* parent, vector<string> label, int nodeCount){
	int root = nodeCount -1;
	int p = parent[node];;
	while(p != root && label[p]==""){
		p = parent[p];
	}
	return p;
}

std::string getGraphVizBinTree(int* parents, int nodeCount, int m, vector<string> label){
	std::stringstream content;
	content << "digraph G {\n";
	content << "node [color=deeppink4, style=filled, fontcolor=white, fontsize=20, fontname=Verdana];\n";
	for(int i=m; i<nodeCount-1; i++){
		if(label[i] != ""){
			int labelledLCA = getLcaWithLabel(i, parents, label, nodeCount);
			content << label[labelledLCA] << " -> " << label[i] << ";\n";
//		if(label[parents[i]] == ""){
//			content  << parents[i] << " -> ";
//		}
//		else{
//			content << label[parents[i]] << " -> ";
//		}
//		if(label[i] == ""){
//		  content  << i << ";\n";
//		}
//		else{
//			content << label[i] << ";\n";

		}
	}
	content << "node [color=lightgrey, style=filled, fontcolor=black];\n";
	for(int i=0; i<m; i++){
		int labelledLCA = getLcaWithLabel(i, parents, label, nodeCount);
		content << label[labelledLCA] << " -> " << "s" << i << ";\n";



//		if(label[parents[i]] == ""){
//			content << parents[i] << " -> ";
//		}
//		else{
//			content << label[parents[i]] << " -> ";
//		}
//
//		content << "s" << i << ";\n";


	}
	content <<  "}\n";
	return content.str();
}



string getMutTreeGraphViz(vector<string> label, int nodeCount, int m, int* parent){
	stringstream nodes;
	stringstream leafedges;
	stringstream edges;
	for(int i=0; i<m; i++){
		if(label.at(i) != ""){
			nodes << "s" << i << "[label=\"s" << i << "\"];\n";                 // si [label="si"];
			nodes        << i << "[label=\"" << label.at(i) << "\"];\n";                 //   i [label="i"];
			leafedges << "s" << i << " -> " << i << ";\n";
			edges <<        i << " -> " << getLcaWithLabel(i, parent, label, nodeCount) << ";\n";
		}
		else{
			nodes << i << "[label=\"s" << i << "\"];\n";
			leafedges << i << " -> " << getLcaWithLabel(i, parent, label, nodeCount) << ";\n";
		}
	}

	stringstream str;

	str << "digraph g{\n";
	str << nodes.str();
	str << "node [color=deeppink4, style=filled, fontcolor=white];	\n";
	str << edges.str();
	str << "node [color=lightgrey, style=filled, fontcolor=black];  \n";
	str << leafedges.str();
	str << "}\n";
	return str.str();
}

/* writes the given string to file */
void writeToFile(string content, string fileName){
	ofstream outfile;
	outfile.open (fileName.c_str());
	outfile << content;
	outfile.close();
}

/* creates the content for the GraphViz file from a parent vector, using numbers as node labels (from 1 to n+1) */
std::string getGraphVizFileContentNumbers(int* parents, int n){
	std::stringstream content;
	content << "digraph G {\n";
	content << "node [color=deeppink4, style=filled, fontcolor=white];\n";
	for(int i=0; i<n; i++){
		content << parents[i]+1  << " -> "  << i+1 << ";\n";      // plus 1 to start gene labeling at 1 (instead of 0)
	}
	content <<  "}\n";
	return content.str();
}


/* creates the content for the GraphViz file from a parent vector, using the gene names as node labels */
std::string getGraphVizFileContentNames(int* parents, int n, vector<string> geneNames, bool attachSamples, bool** ancMatrix, int m, double** logScores, int** dataMatrix){
	std::stringstream content;
	content << "digraph G {\n";
	content << "node [color=deeppink4, style=filled, fontcolor=white];\n";

	for(int i=0; i<n; i++){
		content << geneNames[parents[i]] << " -> "  << geneNames[i]  << ";\n";
	}

	if(attachSamples==true){

		content << "node [color=lightgrey, style=filled, fontcolor=black];\n";
		std::string attachment = getBestAttachmentString(ancMatrix, n, m, logScores, dataMatrix, geneNames);
		content << attachment;
	}
	content <<  "}\n";
	return content.str();
}

/* creates the attachment string for the samples, the optimal attachment points are recomputed from scratch based on error log Scores */
std::string getBestAttachmentString(bool ** ancMatrix, int n, int m, double** logScores, int** dataMatrix, vector<string> geneNames){
	bool** matrix = attachmentPoints(ancMatrix, n, m, logScores, dataMatrix);
	std::stringstream a;
	for(int i=0; i<=n; i++){
		for(int j=0; j<m; j++){
			if(matrix[i][j]==true){
				a << geneNames[i] << " -> s" << j << ";\n";
			}
		}
	}
	return a.str();
}

/* This is a re-computation of the best attachment points of the samples to a tree for printing the tree with attachment points */
/*   gets an ancestor matrix and returns a bit matrix indicating the best attachment points of each sample based on the error log scores */
bool** attachmentPoints(bool ** ancMatrix, int n, int m, double** logScores, int** dataMatrix){

    double treeScore = 0.0;
    bool ** attachment = init_boolMatrix(n+1, m, false);
  	for(int sample=0; sample<m; sample++){       // foreach sample
  		double bestAttachmentScore = 0.0;     // currently best score for attaching sample
  		for(int gene=0; gene<n; gene++){   // start with attaching node to root (no genes mutated)
  			bestAttachmentScore += logScores[dataMatrix[sample][gene]][0];
  		}
  		for(int parent=0; parent<n; parent++){      // try all attachment points (genes)
  		    double attachmentScore=0.0;
  		    for(int gene=0; gene<n; gene++){     // sum up scores for each gene, score for zero if gene is not an ancestor of parent, score for one else wise
  		    	attachmentScore += logScores[dataMatrix[sample][gene]][ancMatrix[gene][parent]];
  		    }
  		    if(attachmentScore > bestAttachmentScore){
  		        bestAttachmentScore = attachmentScore;
  		    }
  		}
  		for(int parent=0; parent<n; parent++){      // try all attachment points (genes)
  		 	double attachmentScore=0.0;
  		 	for(int gene=0; gene<n; gene++){     // sum up scores for each gene, score for zero if gene is not an ancestor of parent, score for one else wise
  		 		attachmentScore += logScores[dataMatrix[sample][gene]][ancMatrix[gene][parent]];
  		 	}
  		  	if(attachmentScore == bestAttachmentScore){
  		  		attachment[parent][sample] = true;
  		  	}
  		}
  		bool rootAttachment = true;
  		for(int parent=0; parent<n; parent++){
  			if(attachment[parent][sample] == true){
  				rootAttachment = false;
  				break;
  			}
  		}
  		if(rootAttachment == true){
  			attachment[n][sample] = true;
  		}
  		treeScore += bestAttachmentScore;
  	}
  	return attachment;
}


/* prints all trees in list of optimal trees to the console, first as parent vector, then as GraphViz file */
void printParentVectors(vector<bool**> optimalTrees, int n, int m, double** logScores, int** dataMatrix){
	for(int i=0; i<optimalTrees.size(); i++){
		int* parents = ancMatrixToParVector(optimalTrees[i], n);
		print_intArray(parents,n);
		//print_boolMatrix(attachmentPoints(optimalTrees[i], n, m, logScores, dataMatrix), n, m);
		printGraphVizFile(parents, n);
	}
}


/* prints the GraphViz file for a tree to the console */
void printGraphVizFile(int* parents, int n){
	cout << "digraph G {\n";
	cout << "node [color=deeppink4, style=filled, fontcolor=white];\n";
	for(int i=0; i<n; i++){
		cout << parents[i] << " -> " << i << "\n";
	}
	cout << "}\n";
}

void printSampleTrees(vector<int*> list, int n, string fileName){
	if(list.size()==0){ return;}
	std::stringstream a;
	for(int i=0; i<list.size(); i++){
		for(int j=0; j<n; j++){
			a << list[i][j];
			if(j<n-1){
				a  << " ";
			}
		}
		a << "\n";
	}
	writeToFile(a.str(), fileName);
	cout << "Trees written to: " << fileName;
}

/* prints the score of the tree predicted by the Kim&Simon approach for the given error log scores */
void printScoreKimSimonTree(int n, int m, double** logScores, int** dataMatrix, char scoreType){
	int parent[] = {2, 4, 17, 2, 9, 9, 2, 2, 4, 18, 2, 1, 2, 2, 9, 2, 2, 11};
	double KimSimonScore = scoreTree(n, m, logScores, dataMatrix, scoreType, parent, -DBL_MAX);
	cout.precision(20);
	cout << "KimSimonScore: " << KimSimonScore << "\n";
}


