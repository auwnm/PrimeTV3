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

*/

#include "Layoutrees.h"

#include "../tree/Treeextended.h"
#include "../Parameters.h"
#include "../tree/Node.hh"

LayoutTrees::LayoutTrees(TreeExtended *r,TreeExtended *g,
                         Parameters *p,const GammaMapEx<Node> *gm,
                         const LambdaMapEx<Node> *l)
    :species(r),gene(g),parameters(p),gamma(gm),lambda(l),
      bv(r->getNumberOfNodes()),Adress(g->getNumberOfNodes())
{


}

void LayoutTrees::start()
{
    if (parameters->ladd == 'r')
    {
        Ladderize_right();
        equal = parameters->equalTimes;
    }
    else if (parameters->ladd == 'l')
    {
        Ladderize_left();
        equal = false;
    }
    
    parameters->maxLeafNameSize = biggestLabel() * parameters->fontsize;
    // we calculate the separation between the tree and the margin of the canvas
    // according to the size of the picture and the size of the biggest leaf
    parameters->separation = (parameters->width / 10) + parameters->maxLeafNameSize;
    parameters->root_sep = (parameters->separation / 2) - parameters->maxLeafNameSize;
    //scale by time and do it equally distributed
    nodetime = parameters->scaleByTime;
    
    calculateSizes();
    calculateIntervals();

    currentY = YCanvasSize - (yspace / 2.0);
    
    CountSpeciesCoordinates(species->getRootNode(),0);
    FindDuplications(gene->getRootNode());
    MapDuplications(gene->getRootNode(),species->getRootNode()->getNumber()+1);
    
    
    initPrePlaces();
    getPlace(species->getRootNode(),gene->getRootNode()); 
    std::cout << "finshed getPlace\n";
    downwardPlaces(gene->getRootNode(),species->getRootNode(), -1);
    std::cout << "finshed downwardPlaces\n";
    setPlaces();
    std::cout << "finshed setPlaces\n";
    
    //setLossPositions();
    
    /*
    Node *u = gene->getNode(61);
    Node *v = u->getLeftChild();
    std::cout<< " ---- > gene Node " << v->getNumber() << endl;

    Node *lgp = gamma->getLowestGammaPath( *v ) ;
    std::cout<< " ---- > lgp Node " << lgp->getNumber() << endl;  
   
    Node *x = species->getNode(12);
    Node *y = x->getDominatingChild( lgp );
    std::cout<< " ---- > Species Node " << y->getNumber() << endl;
    */



    CountGeneCoordinates(gene->getRootNode());
}


LayoutTrees::~LayoutTrees()
{

}

void LayoutTrees::calculateSizes()
{
    //extra espace from the root to the margin to draw gene nodes mapped to the root
    xCanvasXtra = parameters->root_sep * gamma->getSize(species->getRootNode());
    
    if(parameters->horiz)
    {
        yspace =  (parameters->width - parameters->separation) / species->getNumberOfLeaves() ;
    }
    else
    {
        yspace =  (parameters->height - parameters->separation) / species->getNumberOfLeaves() ;
    }
    
       

    NodeHeight = yspace / 2.5;

    //NodeHeight = yspace / 2.5;

    //max numbers of gene nodes mapped to a Species node
    int maxnodesmapped = MostGenes();

    //check if the tree is too big to fit in the actual canvas size
    if(NodeHeight <  (parameters->min_node_height * maxnodesmapped))
    {
        double yrate = ((parameters->min_node_height * maxnodesmapped) - NodeHeight) * species->getNumberOfLeaves();
        parameters->height += yrate;
        yspace =  (parameters->height - parameters->separation) / species->getNumberOfLeaves();
        NodeHeight = yspace / 2.5;
        //assuming the x ampliation will be equal to the y ampliation
        parameters->width += yrate;
    }
    
    if(parameters->horiz)
    {
        YCanvasSize = parameters->width - parameters->separation;
        XCanvasSize = parameters->height - xCanvasXtra - parameters->separation;
    }
    else
    {
        YCanvasSize = parameters->height - parameters->separation;
        XCanvasSize = parameters->width - xCanvasXtra - parameters->separation;
    }

}


void LayoutTrees::calculateIntervals()
{
    maxdeepleaf = species->getRootNode()->getMaxPathToLeaf() + 1;
    
    if(equal)
    {
        maxdeepleaftimes = maptimes();

        for(int i=0; i<maxdeepleaftimes;++i)
        {
            double ratio = (double)i / (maxdeepleaftimes);
            double value = ratio * XCanvasSize;
            double time = maptime[i];
            numXPositionsTimes.insert(pair<double,double>(time,value));
        }
    }
    
    for(int i=0; i<maxdeepleaf;++i)
    {
        double ratio = (double)i / (maxdeepleaf-1);
        double value = ratio * XCanvasSize;
        numXPositions.push_back(value);
    }

}


static bool sort_double(double u, double v)
{
    return u > v;
}

unsigned LayoutTrees::maptimes()
{
    for (unsigned i = 0; i < species->getNumberOfNodes(); i++)
    {
        Node *n = species->getNode(i);
        double time = n->getNodeTime();
        std::vector<double>::iterator it = std::find(maptime.begin(),maptime.end(),time);
        if(it == maptime.end())
        {
            maptime.push_back(n->getNodeTime());
        }
    }
    std::sort(maptime.begin(),maptime.end(),sort_double);
    return maptime.size();
}

double LayoutTrees::getRightMostCoordinate (Node* o)
{
    if (o->isLeaf())
    {
        return o->getY();
    }
    
    else
    {
        return getRightMostCoordinate(o->getRightChild());
    }
}

double LayoutTrees::getLeftMostCoordinate (Node* o)
{
    if(o->isLeaf())
    {
        return o->getY();
    }
    
    else
    {
        return getLeftMostCoordinate(o->getRightChild());
    }
}


/* yspace has been calculated according to the number of leaves and the height
 * so the cordinate y of each node will be increased by y on the leaves, the x
 * cordinate on the leaves is always the same, the y cordinates of intern
 * nodes is calculated in the midpoint between the right most y and the
 * left most y, the x position is calcuted used the time mapped vectors
 */
void LayoutTrees::CountSpeciesCoordinates(Node *n, int depth)
{
     
    //printf("node %d at depth = %d \n",n->getNumber(),depth);
    if (n->isLeaf())
    {
        n->setY(currentY);
        currentY -= yspace;
        n->setX(XCanvasSize + xCanvasXtra);
    }
    else
    {
        Node *left = n->getLeftChild();
        Node *right = n->getRightChild();

        CountSpeciesCoordinates(left, depth + 1);
        CountSpeciesCoordinates(right, depth + 1);

        double sumyleft = getRightMostCoordinate(left);
        double sumyright = getLeftMostCoordinate(right);
        double yposition = (sumyleft + sumyright) / 2;
        double time = n->getNodeTime();
        double xposition;

        if(nodetime && !equal)
        {
            xposition = ((1-time) * XCanvasSize) + xCanvasXtra;
        }
        else if(nodetime && equal)
        {
            xposition = numXPositionsTimes[time] + xCanvasXtra;
        }
        else
        {
            xposition = numXPositions.at(depth) + xCanvasXtra;
        }

        n->setY(yposition);
        n->setX(xposition);

        CalcLegIntersection(left,right,n);
    }
    
}


int LayoutTrees::MostGenes()
{
    int currentMax = 0;
    for(Node *n = species->getPostOderBegin(); n != NULL; n = species->postorder_next(n))
    {
        int size = gamma->getSize(n);
        if( size > currentMax )
        {
            currentMax = size;
        }
    }
    return currentMax;
}

/////////////////////////////////////////////////////////////////////////


/*
double LayoutTrees::getPosition (Node* n)
{
 	   
    if ( n->isLeaf() ){
        return n->getPosition();
    }
    else if( !(n->isLeaf()) && gamma->isSpeciation(*n)) { 
       Node *left  = n->getLeftChild();
       Node *right = n->getRightChild();
           
         
       Node *spnLeft  = gamma->getLowestGammaPath(*left);
       Node *spnRight = gamma->getLowestGammaPath(*right);

       double lx = spnLeft->getX();
       double rx = spnRight->getX();
                 
       if(lx == rx)
           return ( getPosition(left) + getPosition(right) ) / 2.0 ;
       else if( lx < rx )
           return getPosition(left);
        else
           return getPosition(right);  
      
      }
    printf("OK");
}
*/


int LayoutTrees::getPosition (Node* n)
{
    if ( n->isLeaf() || gamma->isSpeciation(*n) )
    {
        return n->getPosition();
    }
    else
    {
	 Node *left  = n->getLeftChild();
         Node *right = n->getRightChild();
	  
         int leftPos  =  getPosition(left);
         int rightPos =  getPosition(right);

         if(leftPos < rightPos)
           return leftPos;
         else
           return rightPos; 	
    }
}


int LayoutTrees::searchPosition(int arr[],int size,int nodeid) {

   for(int i=0; i<size ; i++) {
      if(arr[i] == nodeid) 
	return i;
   }    
   return -1; 

}


int LayoutTrees::getFixedPosition(Node *n) {
   Node *spn = gamma->getLowestGammaPath(*n);
   int size = gamma->getSize(spn);
   int snode = spn->getNumber(); 


    /*	
    switch(snode) {

    case 2:
         return searchPosition(pos2,sizeof(pos2),n->getNumber());
         break;
    
    case 4:
         return searchPosition(pos4,sizeof(pos4),n->getNumber());
         break;
      
    case 6:
         return searchPosition(pos6,sizeof(pos6),n->getNumber());
	 break;

    default:
          return -1;
          break; 
   } 

   */
    return -1;
}


//////////////////////////////////////////////////////////////////////////

void LayoutTrees::initPrePlaces() {
   int size = gene->getNumberOfNodes();
   for(unsigned int u = 0; u < species->getNumberOfNodes(); u++) {
     Node* spn = species->getNode(u);
     spn->initializePrePlaces(size);
     spn->initializePlaces(gamma->getSize(spn));	
  }
}

// Node x is a speciation node
// Node u is a gene node  
std::vector<int> LayoutTrees::getPlace (Node *x, Node *u)
{
  SetOfNodesEx<Node> gamma_set = gamma->getFullGamma(*x);
  if(gamma_set.member(u)) // u is mapping to this edge
    {
      // std::cout<< "host node " << x->getNumber() << " ";
      if(u->isLeaf() && x->isLeaf() && gamma->isInGamma(u, x) )  	// If u is a leaf
	{
	  std::cout<< "host node " << x->getNumber() << ":  guest node "<< u->getNumber() <<" is a leaf at position "<< u->getPosition();
	  std::vector<int> ret(2, -1);
	  ret[0] = u->getPosition();
	  x->prePlace[u->getNumber()] = ret;
	  std::cout << " returning " << ret[0] <<","<< ret[1] <<endl;
	  return ret;
	}
      else if( gamma->isInGamma(u, x) && gamma->isSpeciation(*u) && (*lambda)[u] == x )  // If u is a speciation node
	{
	  std::cout<< "host node " << x->getNumber() << ":  guest node " << u->getNumber() << " is a speciation node, \n";
	  std::vector<int> lpos = getPlace ( x->getLeftChild() , u->getLeftChild()  );
	  std::vector<int> rpos(2,-1);
	  if(lpos != std::vector<int>(2, -1))
	    {
	      std::cout << "host " << x->getNumber() << " guest " << u->getNumber() << "lpos = " << lpos[0] << std::endl;
	      rpos = getPlace ( x->getRightChild() , u->getRightChild() );
	      std::cout << "host " << x->getNumber() << " guest " << u->getNumber() << "rpos = " << rpos[0] << std::endl;
	    }
	  else
	    {
	      lpos = getPlace ( x->getRightChild() , u->getLeftChild()  );
	      rpos = getPlace ( x->getLeftChild() , u->getRightChild() );
	    }
	  if(lpos == std::vector<int>(2, -1) || rpos == std::vector<int>(2, -1)) // check sanity -- This should be an assert in final version
	    {
	      ostringstream oss;
	      oss << "Programming Error: getPlace(" << x->getNumber()
		  << "," << u->getNumber() << "). "
		  << "Inconsistency for child assumptions in speciaiton case";
	      throw AnError(oss.str());
	    }
	  
	  std::vector<int> ret = intersectChildPlaces(lpos,rpos);
	  x->prePlace[u->getNumber()] = ret;
	  std::cout<< "host node " << x->getNumber() << ":  guest node " << u->getNumber() << " is a speciation node, returning " << ret[0] <<"," << ret[1] << endl;	
	  return ret; 
	}
      else if(gamma_set.member(u->getLeftChild()) && gamma_set.member(u->getRightChild())) // u is a duplication on edge to x
	{  
	  
	  std::cout<< "host node " << x->getNumber() << ":  guest node " << u->getNumber() << " is a duplication node,\n";
	  std::vector<int> lpos = getPlace ( x, u->getLeftChild()  );
	  std::vector<int> rpos = getPlace ( x, u->getRightChild() );
	  std::vector<int> ret = intersectChildPlaces(lpos,rpos);
	  
	  std::cout<< "host node " << x->getNumber() << ":  guest node " << u->getNumber() << " is a duplication node, returning " << ret[0]<<","<<ret[1] << endl;
	  
	  return ret; 	  
	}
      else    // for loss node
	{
	  std::cout<< "host node " << x->getNumber() << ":  guest node " << u->getNumber() << " is a loss, \n";
	  std::vector<int> pos(2,-1);
	  // Either lpos or rpos should return -1 since u can't map to both
	  std::vector<int> lpos = getPlace (x->getLeftChild(), u);
	  std::vector<int> rpos = getPlace (x->getRightChild(), u);
	  if(lpos == std::vector<int>(2, -1))
	    {
	      if(rpos == std::vector<int>(2, -1)) // check sanity -- this should be an assert in final version
		{ 
		  ostringstream oss;
		  oss << "Programming Error: getPlace(" << x->getNumber() << "," << u->getNumber() << ")";
		  throw AnError(oss.str());
		}
	      pos = rpos;
	    }
	  else
	    {
	      if(rpos != std::vector<int>(2, -1)) // check sanity -- this should be an assert in final version
		{ 
		  ostringstream oss;
		  oss << "Programming Error: getPlace(" << x->getNumber() << "," << u->getNumber() << ")";
		  throw AnError(oss.str());
		}
	      pos = lpos;
	    }
	  
	  x->prePlace[u->getNumber()]	= pos;
	  std::cout<< "host node " << x->getNumber() << ":  guest node " << u->getNumber() << " is a loss, returning "<< pos[0] << ","<<pos[1] << endl;

	  return pos;
	}  
    }
  else
    {
      return std::vector<int>(2, -1); // u is not mapped to the subtree rooted at x
    }
  
}
// Node x is a speciation node
// Node u is a gene node  
// int LayoutTrees::getPlace (Node *x, Node *u)
// {
//   SetOfNodesEx<Node> gamma_set = gamma->getFullGamma(*x);
//   if(gamma_set.member(u)) // u is mapping to this edge
//     {
//       std::cout<< "host node " << x->getNumber() << " ";
//       if(u->isLeaf() && x->isLeaf() && gamma->isInGamma(u, x) )  	// If u is a leaf
// 	{
// 	  std::cout<< "guest node "<< u->getNumber() <<" is a leaf at position "<< u->getPosition() << endl;
// 	  x->prePlace[u->getNumber()] = u->getPosition();
// 	  return u->getPosition();	
// 	}
//       else if( gamma->isInGamma(u, x) && gamma->isSpeciation(*u) && (*lambda)[u] == x )  // If u is a speciation node
// 	{
// 	  std::cout<< "guest node " << u->getNumber() << " is a speciation node." << endl;	
	  
// 	  int lpos = getPlace ( x->getLeftChild() , u->getLeftChild()  );
// 	  int rpos = -1;
// 	  if(lpos != -1)
// 	    {
// 	      rpos = getPlace ( x->getRightChild() , u->getRightChild() );
// 	    }
// 	  else
// 	    {
// 	      lpos = getPlace ( x->getRightChild() , u->getLeftChild()  );
// 	      rpos = getPlace ( x->getLeftChild() , u->getRightChild() );
// 	    }
// 	  if(lpos == -1 || rpos == -1) // check sanity -- This should be an assert in final version
// 	    {
// 	      ostringstream oss;
// 	      oss << "Programming Error: getPlace(" << x->getNumber()
// 	  	  << "," << u->getNumber() << "). "
// 	  	  << "Inconsistency for child assumptions in speciaiton case";
// 	  	throw AnError(oss.str());
// 	    }
	  
// 	  std::cout<< "guest node " << u->getNumber() << " lpos = " << lpos << " rpos = " << rpos << "\n";
	  
// 	  x->prePlace[u->getNumber()] = std::min(lpos,rpos);
// 	  return std::min(lpos,rpos);  	
// 	}
//       else if(gamma_set.member(u->getLeftChild()) && gamma_set.member(u->getRightChild())) // u is a duplication on edge to x
// 	{  
// 	  std::cout<< "guest node " << u->getNumber() << " is a duplication node." << endl;
// 	  int lpos = getPlace ( x, u->getLeftChild()  );
// 	  int rpos = getPlace ( x, u->getRightChild() );
// 	  return std::min(lpos,rpos);
	  
// 	}
//       else    // for loss node
// 	{
// 	  std::cout<< "guest node " << u->getNumber() << " is a loss." << endl;
	  
// 	  int pos = -1;
// 	  // Either lpos or rpos should return -1 since u can't map to both
// 	  int lpos = getPlace (x->getLeftChild(), u);
// 	  int rpos = getPlace (x->getRightChild(), u);
// 	  if(lpos == -1)
// 	    {
// 	      if(rpos == -1) // check sanity -- this should be an assert in final version
// 	  	{ 
// 	  	  ostringstream oss;
// 	  	  oss << "Programming Error: getPlace(" << x->getNumber() << "," << u->getNumber() << ")";
// 	  	  throw AnError(oss.str());
// 	  	}
// 	      pos = rpos;
// 	    }
// 	  else
// 	    {
// 	      if(rpos != -1) // check sanity -- this should be an assert in final version
// 	  	{ 
// 	  	  ostringstream oss;
// 	  	  oss << "Programming Error: getPlace(" << x->getNumber() << "," << u->getNumber() << ")";
// 	  	  throw AnError(oss.str());
// 	  	}
// 	      pos = lpos;
// 	    }	   
	  	  
// 	  x->prePlace[u->getNumber()]	= pos;
// 	  return pos;
// 	}  
//     }
//   else
//     {
//       return -1; // u is not mapped to the subtree rooted at x
//     }
  
// }

std::vector<int>
LayoutTrees::intersectChildPlaces(std::vector<int> lpos, std::vector<int> rpos)
{
  if(lpos[0] != -1 || rpos[0]!=-1) // Sanity
    {
      
      int l = -1;
      int r = -1;
      int min = 999999999;
      for(std::vector<int>::const_iterator i = lpos.begin(); i != lpos.end(); i++)
	{
	  if(*i != -1)
	    {	
	      for(std::vector<int>::const_iterator j = rpos.begin(); j != rpos.end(); j++)
		{
		  if(*j != -1)
		    {
		      int diff = std::abs(*i-*j);
		      if(diff < min)
			{
			  min  = diff;
			  l = *i;
			  r = *j;
			  std::cout << "l = " << l << " r = " << r << " min = " << min << "\n";
			}
		    }
		}
	    }
	}
      if(l == -1)
	{
	  throw AnError("intersectChildPlaces is called with empty pairs");
	}
      std::vector<int> ret;
      ret.push_back(l);
      ret.push_back(r);
      return ret;
    }
  else
    {
      return std::vector<int>(2,-1);
    }
}

void
LayoutTrees::downwardPlaces(Node* u, Node* x, int prev)
{
  std::cout << "downwardPlaces(" << u->getNumber() <<"," <<x->getNumber() <<"," << prev<< ")" << std::endl;
  SetOfNodesEx<Node> gamma_set = gamma->getFullGamma(*x);
  if(gamma_set.member(u)) // u is mapping to this edge
    {
      std::cout  << "downwardPlaces(" << u->getNumber() <<"," <<x->getNumber() <<"," << prev<< ")"<< " is in gamma" << std::endl;
      std::vector<int>& tmp = x->prePlace[u-> getNumber()];
      if(tmp[0] != -1)
	{
	  if(u->isLeaf() && x->isLeaf() && gamma->isInGamma(u, x) )  	// If u is a leaf
	    {
	      std::cout << "downwardPlaces(" << u->getNumber() <<"," <<x->getNumber() <<"," << prev<< ")"<< " a leaf, untoouched " << tmp[0] << std::endl;
	      return;
	    }
	  else 
	    {
	      std::cout << "downwardPlaces(" << u->getNumber() <<"," <<x->getNumber() <<"," << prev<< ")" << " is in psec or loss"<< std::endl;
	      if(prev != -1)
		{
		  if(std::abs(prev-tmp[1]) < std::abs(prev-tmp[0]))
		    {
		      prev = tmp[1];
		      std::reverse(tmp.begin(), tmp.end());
		      x->prePlace[u-> getNumber()] = tmp;
		      std::cout << "downwardPlaces(" << u->getNumber() <<"," <<x->getNumber() <<"," << prev<< ")"<< "tmp0 = " << tmp[0] << " tmp1 = " <<tmp[1] << " change to " << tmp[0] << std::endl;
		    }
		  else
		    {
		      std::cout << "downwardPlaces(" << u->getNumber() <<"," <<x->getNumber() <<"," << prev<< ")" << " untouched " << tmp[0] << std::endl;
		    }
		}
	      else
		{
		  prev = tmp[0];
		  std::cout  << "downwardPlaces(" << u->getNumber() <<"," <<x->getNumber() <<"," << prev<< ")"<< " is on host root, untouched " << tmp[0] << std::endl;
		}
	      if( gamma->isInGamma(u, x) && gamma->isSpeciation(*u) && (*lambda)[u] == x )  // If u is a speciation node
		{
		  downwardPlaces(u->getLeftChild(), x->getLeftChild(), prev);
		  downwardPlaces(u->getRightChild(), x->getRightChild(), prev);
		  downwardPlaces(u->getRightChild(), x->getLeftChild(), prev);
		  downwardPlaces(u->getLeftChild(), x->getRightChild(), prev);
		}
	      else //loss
		{
		  downwardPlaces(u, x->getLeftChild(), prev);
		  downwardPlaces(u, x->getRightChild(), prev);
		}
	    }
	}
      else
	{
	  std::cout  << "downwardPlaces(" << u->getNumber() <<"," <<x->getNumber() <<"," << prev<< ")"<< " is duplications " << std::endl;

	  downwardPlaces(u->getLeftChild(), x, prev);
	  downwardPlaces(u->getRightChild(), x, prev);
	  return;
	}
    }
}
  
void LayoutTrees::setPlaces()
{
  int size = gene->getNumberOfNodes();
  std::vector < nodeplace > nvec;
 
  bool showPlaces = true;

  for(unsigned int x = 0; x < species->getNumberOfNodes(); x++) {
     Node* spn = species->getNode(x);
     int sizeX = gamma->getSize(spn);
     
     nvec.clear();
     if(showPlaces) std::cout<< " Species Node " << spn->getNumber() << " : [ ";
     for (unsigned i = 0; i < gene->getNumberOfNodes(); i++) {
       int pos = spn->prePlace[i][0];
        // if( spn->prePlace[i] == -1 ) continue;
       if( pos == -1 ) continue;
       // if(showPlaces) std::cout<< "( "<< i << ", "<< spn->prePlace[i] << ") ";  
       if(showPlaces) std::cout<< "setPlaces()" << " host node " << x <<"( "<< i << ", "<< pos << ") " << std::endl;  
        // nvec.push_back( nodeplace(i,  double(spn->prePlace[i])/sizeX) );
        nvec.push_back( nodeplace(i,  double(pos)/sizeX) );
     }
     if(showPlaces) std::cout<< " ]" << endl;


     // std::sort(nvec.begin(), nvec.end());
     std::stable_sort(nvec.begin(), nvec.end());
     for(std::vector<nodeplace>::size_type i = 0; i != nvec.size(); i++) {
         spn -> setPlace(i,nvec[i].id);
     }
    
  }

  /* 
  Node *u = gene->getNode(84);
  Node *top  = gamma->getHighestGammaPath(*u);
  // Node *down = gamma->getLowestGammaPath(*u);
  std::cout<< " ---> LGP Node " << top->getNumber() << endl;
  */


}

void LayoutTrees::AssignLeafGene(Node *n)
{
    Node *spn = gamma->getLowestGammaPath(*n);
    n->setX(spn->getX());
    double y;
    int size = gamma->getSize(spn);

    if(size > 1)
    {
       int yoffset = spn->searchPlace(n->getNumber());

       if(yoffset == -1)
       	  yoffset = spn->getVisited();

	int delta = ( 2 * NodeHeight / (size - 1) );    
 	y = (spn->getY() + NodeHeight) - (delta * yoffset);
		
    }
    else
    {
        y = spn->getY();
    }
    
    spn->incVisited();
    n->setY(y);
    n->setHostChild(spn);
    
    if(!n->isRoot() && gamma->isLateralTransfer(*n->getParent())
            && ((*lambda)[n] == (*lambda)[n->getParent()]))
    {
        Node *destiny = (*lambda)[getHighestMappedLGT(n)];
        n->setHostParent(destiny);
    }
    else
    {
        n->setHostParent(spn);
    }
}


/////////////////////////////////////////////////////////////////////////
/* each species node has a visited atribute, if the size (number os nodes mapped)
 * of a species node is > 1 then the nodeheight will be divided by the number
 * of nodes mapped and the current node located in the according position
 */
/*
void LayoutTrees::AssignLeafGene(Node *n)
{
    Node *spn = gamma->getLowestGammaPath(*n);
    n->setX(spn->getX());
    double y;
    int size = gamma->getSize(spn);
    
    if(size > 1)
    {
        int yoffset = spn->getVisited();         // every gene node has its own offset
        int delta = ( NodeHeight / (size - 1) ); // +2 ;    // Changed delta is space between node
        y = (spn->getY() - NodeHeight/2) + (delta * yoffset);
	
    }
    else
    {
        y = spn->getY();
    }
    
    spn->incVisited();
    n->setY(y);
    n->setHostChild(spn);
    
    if(!n->isRoot() && gamma->isLateralTransfer(*n->getParent())
            && ((*lambda)[n] == (*lambda)[n->getParent()]))
    {
        Node *destiny = (*lambda)[getHighestMappedLGT(n)];
        n->setHostParent(destiny);
    }
    else
    {
        n->setHostParent(spn);
    }
}
*/

void LayoutTrees::CountGeneCoordinates(Node* n)
{
    if(n->isLeaf())
    {
        n->setReconcilation(Leaf);
        AssignLeafGene(n);
	
    }
    else
    {
        Node *left = n->getLeftChild();
        Node *right = n->getRightChild();
        
        CountGeneCoordinates(left);
        CountGeneCoordinates(right);

        if(gamma->isSpeciation(*n) && !gamma->isLateralTransfer(*n)) //speciation
        {
            n->setReconcilation(Speciation);
            AssignLeafGene(n);
        }
        else if (gamma->isLateralTransfer(*n)) //lateral transfer
        {
            AssignGeneLGT(n);
        }
        else //duplication
        {
            AssignGeneDuplication(n);
        }
    }
}

void LayoutTrees::AssignGeneDuplication(Node *n)
{
  Node *spb = Adress[n]; 
  double proportion = 0; 
  double delta = 0; 
  double edge = 0; 
  if(!spb->isRoot())
    {
      Node *spbP = spb->getParent();
      edge = spb->getX() - spbP->getX();
      proportion = ((spbP->getY() - spb->getY()) / edge);
      n->setHostParent(spbP);
    }
  else
    {
      edge = spb->getX();
    }
  
  /* we obtain the number of duplication and the duplications levels
   * to figure out the x position of the node, we use the left most and right
   * most cordinates of the duplication to figure out the y position
   */
  double ndupli = bv[spb]+1;
  unsigned duplilevel = Duplevel(n,spb->getNumber()); // I think this is the number of dups below and including n on spb
  delta = (edge/ndupli)*duplilevel;

  cout << "Node " << n->getNumber() << ": ndupli = " << ndupli << ",  duplilevel = " << duplilevel << ", delta = " << delta << "and proportion = " << proportion << "\n";
  
  n->setX(spb->getX()-delta);
  
  
  double rightMost = RightMostCoordinate(n,spb,duplilevel);
  double leftMost  = LeftMostCoordinate(n,spb,duplilevel);

  if(spb->isRoot())
    {
      n->setY( leftMost +  (proportion * delta) );
    }
  else
    {
      double leftChildPlace = getRightMostChildPlace(n->getLeftChild(),spb);
      double rightChildPlace = getLeftMostChildPlace(n->getRightChild(),spb);
      double parentPlace = getParentPlace(n, spb->getParent());
      
      
      if(abs(rightChildPlace - parentPlace) < abs(leftChildPlace - parentPlace))
	{
	  cout << "AssignGeneDuplication(guest Node " << n->getNumber() << ") with host node "<< spb->getNumber() << ": leftChildPlace = " << leftChildPlace << ",  rightChildPlace = " << rightChildPlace << ", parentPlace(" << n->getNumber() <<","<<spb->getParent()->getNumber() << " = " << parentPlace << " chose rightChild" <<std::endl;
	  n->setY( rightMost +  (proportion * delta) );
	}
      else
	{
	  cout << "AssignGeneDuplication(guest Node " << n->getNumber() << ") with host node "<< spb->getNumber() << ": leftChildPlace = " << leftChildPlace << ",  rightChildPlace = " << rightChildPlace << ", parentPlace(" << n->getNumber() <<","<<spb->getParent()->getNumber() << " = " << parentPlace << " chose leftChild" <<std::endl;
	  n->setY( leftMost +  (proportion * delta) ); 
	} 
    }    
  // n->setY( ((rightMost + leftMost) /2) +  (proportion * delta) );
  
  // n->setY( rightMost +  (proportion * delta) ); 
  
  
  n->setReconcilation(Duplication);
  n->setHostChild(spb);
}

void LayoutTrees::AssignGeneLGT(Node *n)
{
    n->setReconcilation(LateralTransfer);
    
    Node *SoriginLT = (*lambda)[n];
    Node *SdestinyLT = (*lambda)[n->getLeftChild()];
    
    if(SoriginLT == SdestinyLT)
    {
        SdestinyLT = (*lambda)[n->getRightChild()];
    }
    
    n->setHostParent(SdestinyLT);
    n->setHostChild(SoriginLT);
}

Node*
LayoutTrees::FindDuplications(Node* node) 
{
    if(gamma->isSpeciation(*node) || gamma->isLateralTransfer(*node))
    {
        if(!node->isLeaf())
        {
            FindDuplications(node->getLeftChild());
            FindDuplications(node->getRightChild());
        }
        
        Node *top = gamma->getHighestGammaPath(*node);
        return top;
    }
    else
    {
        Node *top_dup_l = FindDuplications(node->getLeftChild());
        Node *top_dup_r = FindDuplications(node->getRightChild());
        Node *top_l = gamma->getHighestGammaPath(*(node->getLeftChild()));
        Node *top_r = gamma->getHighestGammaPath(*(node->getRightChild()));

        if (top_l == NULL)
        {
            top_l = top_dup_l;
        }
        if (top_r == NULL)
        {
            top_r = top_dup_r;
        }
        
        Adress[node] = species->mostRecentCommonAncestor(top_l, top_r);
        return Adress[node];
    }
}

Node*
LayoutTrees::MapDuplications(Node* de, unsigned line) 
{
    if(de->isLeaf())
    {
        return de->getParent();
    }

    else
    {
        if(!gamma->isSpeciation(*de) && !gamma->isLateralTransfer(*de))
        {
            unsigned testlineage = Adress[de]->getNumber();

            if(testlineage != line)
            {
                Node* nd = Adress[de];
                unsigned n_levels = Duplevel(de,testlineage);
                if (n_levels > bv[nd])
                {
                    bv[nd] = n_levels;
                }
            }

            line = testlineage;
        }
        
        Node* newNode = MapDuplications(de->getLeftChild(), line);
        Node* newNode2 = MapDuplications(newNode->getRightChild(), line); // Is this aimed to trap LGTs or what? How?

        return newNode2->getParent();
    }
}


// It appears that the aim of this function is to count the number of duplications occuring
// below nd on the same edge (which here is identified by its id-number) including the duplication at nd.
unsigned 
LayoutTrees::Duplevel(Node* nd, int levellineage)                                                          
{
  // Here, second test is whether levellineage is lca of nd in host tree
  // However, all calls to Duplevel use Adress[nd]-> getNumber() as levellineage so the
  // second part of the test below is meaningless!!!!
  if(gamma->isSpeciation(*nd) || Adress[nd]->getNumber() != unsigned(levellineage)) 
    {
        return 0;
    }
    else
    {
        int left = Duplevel(nd->getLeftChild(), levellineage);
        int right = Duplevel(nd->getRightChild(), levellineage);
        return max(left, right) + 1;
    }
}

// This shuld return the leftmost y-coordinate of o or any descendant of o at the vertex end-of_slice
double 
LayoutTrees::LeftMostCoordinate(Node* o, Node *end_of_slice, int duplevel)
{
    if (gamma->isSpeciation(*o))
    {
      if(gamma->getLowestGammaPath(*o) != end_of_slice)  
        {
	  // we should use searchPlaces here
	  return getYforLosses(o, end_of_slice);

	  // double size = gamma->getSize(end_of_slice);
	  // double delta = NodeHeight / size - 1;
	  // double y = (end_of_slice->getY() - NodeHeight/2) + ((duplevel) * delta);
	  // return y;
        }
        else
        {
            return o->getY();
        }
    }
    else
    {
        if (end_of_slice == Adress[o])
        {
            return LeftMostCoordinate(o->getLeftChild(), end_of_slice,duplevel);
        }
        else
        {
	  // we should use searchPlaces here
	  return getYforLosses(o, end_of_slice);

	  // double size = gamma->getSize(end_of_slice);
          //   double delta = NodeHeight / size - 1;
          //   double y = (end_of_slice->getY() - NodeHeight/2) + ((duplevel) * delta);
          //   return y;
        }
    }
}


// This shuld return the right most y-coordinate of o or any descendant of o at the vertex end-of_slice
double
LayoutTrees::RightMostCoordinate (Node* o, Node *end_of_slice, int duplevel)
{
    if (gamma->isSpeciation(*o))
    {
      if(gamma->getLowestGammaPath(*o) != end_of_slice) // o loss at end-of-slice
        {
	  // we should use searchPlaces here
	  return getYforLosses(o, end_of_slice);
	  
            // double size = gamma->getSize(end_of_slice);
            // double delta = NodeHeight / size - 1;
            // double y = (end_of_slice->getY()- NodeHeight/2) + ((duplevel) * delta);
            // return y;
        }
        else
        {
            return o->getY();
        }
    }
    else
    {
      if (end_of_slice == Adress[o]) // There must be a node below o in end_of_slice edge
        {
            return RightMostCoordinate(o->getRightChild(), end_of_slice,duplevel);
        }
      else // o loss at end-of-slice
        {
	  // we should use searchPlaces here
	  return getYforLosses(o, end_of_slice);
	  
	  // double size = gamma->getSize(end_of_slice);
	  // double delta = NodeHeight / size - 1;
	  // double y = (end_of_slice->getY()- NodeHeight/2) + ((duplevel) * delta);
	  // return y; ;
        }
    }
}

// This shuld return the y-coordinate of u or the closest left-most child of u mapping at vertex x
double
LayoutTrees::getLeftMostChildPlace(Node* u, Node *x)
{
  int p = x->searchPlace(u->getNumber());
  if(p == -1)
    {
      return double(getLeftMostChildPlace(u->getLeftChild(), x))/gamma->getSize(x);
    }
  else
    {
      return double(p)/gamma->getSize(x);
    }
}

// This shuld return the y-coordinate of u or the closest right-most child of u mapping at vertex x
double
LayoutTrees::getRightMostChildPlace(Node* u, Node *x)
{
  int p = x->searchPlace(u->getNumber());
  if(p == -1)
    {
      return double(getRightMostChildPlace(u->getRightChild(), x))/gamma->getSize(x);
    }
  else
    {
      return double(p)/gamma->getSize(x);
    }
}

// This shuld return the y-coordinate of u or the closest parent of u mapping at vertex x
double
LayoutTrees::getParentPlace(Node* u, Node *x)
{
  int p = x->searchPlace(u->getNumber());
  if(p == -1)
    {
      if(u->getParent()->isRoot())
	{
	  ostringstream oss;
	  oss << "getParentPlace(" << u->getNumber()
	      << "," << x->getNumber()
	      << "): u->getParent() is called and returns "
	      << u->getParent()->getNumber() << ", which is root\n";
	  throw AnError(oss.str());
	}
      std::cout << "getParentPlace(" << u->getNumber() <<";" << x->getNumber() << ") = " <<getParentPlace(u->getParent(), x)  << "/" << gamma->getSize(x)<< "\n";
      return getParentPlace(u->getParent(), x);
    }
  else
    {
      double ret = double(p)/gamma->getSize(x);
      std::cout << "getParentPlace(" << u->getNumber()<< ";" << x->getNumber() << ") = " <<p  << "/" << gamma->getSize(x)<< " = " << ret << "\n";
      return ret; // double(p)/gamma->getSize(x);
    }
}


// This returns the y-ccordinate of a loss u at x
double
LayoutTrees::getYforLosses(Node* u, Node *x)
{
  int y = -1;
  double size = gamma->getSize(x);
  if(size > 1)
    {
      double delta = NodeHeight / size - 1;
      int yoffset = x->searchPlace(u->getNumber());
      
      if(yoffset == -1) //needed?
	yoffset = x->getVisited();
      
      y = (x->getY() + NodeHeight/2) - (delta * yoffset);
    }
  else
    {
      y = x->getY();
    }
  return y;
}


double LayoutTrees::getNodeHeight()
{
    return NodeHeight;
}

int
LayoutTrees::Ladderize_left() 
{
    return Ladderize_left(species->getRootNode());
}

int
LayoutTrees::Ladderize_left(Node* n) 
{
    if(n->isLeaf())
    {
        return 1;
    }
    else
    {
        int leftsize = Ladderize_left(n->getLeftChild());
        int rightsize = Ladderize_left(n->getRightChild());
        if(leftsize > rightsize)
        {
            n->rotate();
        }
        return leftsize + rightsize;
    }
}

int
LayoutTrees::Ladderize_right() 
{
    return Ladderize_right(species->getRootNode());
}

int
LayoutTrees::Ladderize_right(Node* n)
{
    if(n->isLeaf())
    {
        return 1;
    }
    else
    {
        int leftsize = Ladderize_right(n->getLeftChild());
        int rightsize = Ladderize_right(n->getRightChild());
        if(rightsize > leftsize)
        {
            n->rotate();
        }
        return leftsize + rightsize;
    }
}


void
LayoutTrees::CalcLegIntersection(Node *left, Node *right, Node *u)
{
    double x0, y0, x1, y1, x2, y2, x3, y3;

    x0 = left->getX();
    y0 = left->getY() - NodeHeight;

    x1 = u->getX();
    y1 = u->getY() - NodeHeight;

    x2 = right->getX();
    y2 = right->getY() + NodeHeight;

    x3 = u->getX();
    y3 = u->getY() + NodeHeight;
    
    // The slants of the two lines
    double k_L = (y1 - y0) / (x1 - x0);
    double k_R = (y3 - y2) / (x3 - x2);

    double D_R = (y3 - y0 - k_L * (x3 - x0)) / (k_L - k_R);
    double D_L = x3 + D_R - x0;

    u->setX(x0 + D_L);
    u->setY(y0 + k_L * D_L);
}

//get the highest not LGT mapped node of n
Node* LayoutTrees::getHighestMappedLGT(Node *n)
{
    Node *parent = n->getParent();

    while(gamma->isLateralTransfer(*parent) && !parent->isRoot())
        parent = parent->getParent();

    while(!species->descendant((*lambda)[n],(*lambda)[parent]) && !parent->isRoot())
        parent = parent->getParent();

    return parent;
}

double LayoutTrees::biggestLabel()
{
    double size = 0.0;
    
    for(unsigned i = 0; i < gene->getNumberOfNodes(); i++)
    {
        Node *n = gene->getNode(i);
        if(n->isLeaf())
        {
            size = max(size,(double)n->getName().size());
        }
    }
    
    for(unsigned i = 0; i < species->getNumberOfNodes(); i++)
    {
        Node *n = species->getNode(i);
        if(n->isLeaf())
        {
            size = max(size,(double)n->getName().size());
        }
    }
    
    return size;
}

void LayoutTrees::replaceNodes(const std::map<int,int>& replacements)
{
    for(std::map<int,int>::const_iterator it = replacements.begin();
        it != replacements.end(); ++it)
    {
        Node *first = gene->getNode(it->first);
        Node *second = gene->getNode(it->second);
        std::cout << "Replacing Node " << it->first << " with Node " << it->second << std::endl;
        if(first && second)
        {
            //HostParent and HostChild should not change
            double temp_x = first->getX();
            double temp_y = first->getY();
            first->setX(second->getX());
            first->setY(second->getY());
            second->setX(temp_x);
            second->setY(temp_y);
        }
    }
}










