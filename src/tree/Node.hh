/*
    PrimeTV2 : a visualizer for phylogenetic reconciled trees.
    Copyright (C) 2011  <Jose Fernandez Navarro> <jc.fernandez.navarro@gmail.com>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
    
    Author : Jose Fernandez Navarro  -  jc.fernandez.navarro@gmail.com
             Lars Arvestad, © the MCMC-club, SBC, all rights reserved
 */

#ifndef NODE_HH
#define NODE_HH

#include <iostream>
#include <string>
#include <map>
#include <vector>

#include "../draw/Color.h"
#include "../reconcilation/SetOfNodesEx.hh"

using namespace std;
class Tree;

enum  Type{Leaf=0, Speciation=1, Duplication=2, LateralTransfer=3, Undefined=4};
typedef float Real;

class Node
{

public:

    //Auwn variables
    int position;
    int * place;
    // int * prePlace;
  std::vector< std::vector<int> >  prePlace;
    int numberOfPlaces; 	


    Node(unsigned id);
    Node(unsigned id, const std::string& nodeName);
    virtual ~Node();
    Node(const Node &);

    //Extra Methods aded by Jose Carlos Fernandez
    void setColor(Color c);
    void setSize(double s);
    void setX(double x);
    void setY(double y);
    void setHostParent( Node *hostparent);
    void setHostChild( Node *hostchild);
    void setReconcilation(Type t);
    void setVisited(unsigned inc);
    void incVisited();
    Color getColor();
    double getSize();
    double getX();
    double getY();
    double getX() const;
    double getY() const;
    Node *getHostParent();
    Node *getHostChild();
    Type getReconcilation();
    Node *getHostParent() const;
    Node *getHostChild() const;
    Type getReconcilation() const;
    unsigned getVisited();
    void setRightChild(Node *);
    void setLeftChild(Node *);
    //Extra Methods aded by Jose Carlos Fernandez

    void rotate();
    Node* getLeftChild() const;
    Node* getRightChild() const;
    Node* getParent() const;
    Node* getSibling() const;
    Node* getDominatingChild(Node* y);
    const std::string& getName() const;
    Tree* getTree();
    unsigned getNumber() const;
    unsigned getPorder() const;
    unsigned getNumberOfLeaves() const;
    Real getBranchLength() const;
    Node& operator=(const Node &);
    unsigned getMaxPathToLeaf();
    void setName(const std::string& nodeName);
    void setTree(Tree& T);
    void setChildren(Node *left, Node *right);
    void setParent(Node *parent);
    void changeID(unsigned newID);
    void deleteSubtree();
    SetOfNodesEx<Node> getLeaves();
    bool isLeaf() const;
    bool isRoot() const;
    bool operator<=(const Node& b) const;
    bool operator<(const Node& b) const;
    bool operator<(const Node* b) const;
    bool operator>(const Node& b) const;
    bool dominates(const Node &b) const;
    bool strictlyDominates(const Node &b) const;

    Real getNodeTime() const;
    Real getTime() const;
    Real getLength() const;
    bool changeNodeTime(const Real &t);
    bool changeTime(const Real &et);
    void setNodeTime(const Real &t);
    void setTime(const Real &t);
    void setLength(const Real &newLength);
    
    friend std::ostream& operator<< (std::ostream& o, const Node &v);
    friend std::ostream& operator<< (std::ostream& o, const Node *v);

    //Extra Methods aded by Auwn
    int getPosition();
    void setPosition(int position);  		
   
    void initializePlaces(int size); 
    void initializePrePlaces(int size); 

    bool setPlace(int idx,int nodeid);
    int getPlace(int idx);
    int searchPlace(int nodeid);
    int getNumberOfPlaces();
    int * lossPlace;
private:

    std::string stringify(std::string tag, Real val) const;
    std::string stringify(std::string tag, std::string val) const;
    std::string stringify(std::string tag, Node *v) const;

protected:

    unsigned number;      // the number, the unique id in the tree.

    Node *parent;         // pointer to the parent Node
    Node *leftChild;      // pointer to the left child Node
    Node *rightChild;     // pointer to the right child Node
    long porder;          // Defining partial order of tree.

    Real time;		      // the arc time from the parent to current node
    Real nodeTime;        // the time from the leaves.
    Real branchLength;	  // equals time * rate, we might want to remove this

    std::string name;     // the (leaf) name
    Tree* ownerTree;      // The tree to which I belong

    /* extra features Jose Carlos Fernandez */
    Color color;
    double size;
    double x;
    double y;
    Node *hostParent;
    Node *hostChild;
    Type reconcilation;
    unsigned visited;
    

      
};

#endif
