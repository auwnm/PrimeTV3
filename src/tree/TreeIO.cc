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
             Bengt Sennblad, © the MCMC-club, SBC, all rights reserved
 */

#include "TreeIO.hh"

#include <cassert>		// For early bug detection
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <cmath>
#include <algorithm>

using namespace std;

TreeIO::TreeIO(enum TreeSource source, const std::string &s)
    : source(source),
      stringThatWasPreviouslyNamedS(s)
{}

// Update: Since I sometimes want to output trees without reading 
// anything, I will now allow instantiating the empty object.
//--------------------------------------------------------------------
TreeIO::TreeIO()
    : source(readFromStdin),
      stringThatWasPreviouslyNamedS("")
{}

TreeIO::~TreeIO()
{}

TreeIO::TreeIO(const TreeIO &io)
    : source(io.source),
      stringThatWasPreviouslyNamedS(io.stringThatWasPreviouslyNamedS)
{}

TreeIO&
TreeIO::operator=(const TreeIO &io)
{
    if (this != &io)
    {
        stringThatWasPreviouslyNamedS = io.stringThatWasPreviouslyNamedS;
        source = io.source;
    }
    return *this;
}

TreeIO
TreeIO::fromFile(const std::string &f)
{
    if (f=="")
    {
        return TreeIO();
    }
    else
    {
        return TreeIO(readFromFile, f);
    }
}


TreeIO
TreeIO::fromString(const std::string &s)
{
    return TreeIO(readFromString, s);
}


// Change source using these utilities:
//--------------------------------------------------------------------
void
TreeIO::setSourceFile(const std::string &filename)
{
    source = readFromFile;
    stringThatWasPreviouslyNamedS = filename;
}

void
TreeIO::setSourceString(const std::string &str)
{
    source = readFromString;
    stringThatWasPreviouslyNamedS = str;
}

//----------------------------------------------------------------------
// Read associations
//----------------------------------------------------------------------

// Map leaves in the gene tree to leaves in the species tree
// This is a bit incongruent with the rest of the code and should 
// probably get its own class! /arve
// Expected line format: 
// <whitespace>? <gene name> <whitespace> <species name> <whitespace>?
// Line length is at most 1024 chars.
//--------------------------------------------------------------------
StrStrMap
TreeIO::readGeneSpeciesInfo(const std::string &filename)
{
    ifstream is(filename.c_str());
    int lineno = 1;

    StrStrMap gene2species;
    if(is.peek() == '#') // gs may start with a '#'
    {
        char dummy[LINELENGTH];
        is.getline(dummy, LINELENGTH);// >> dummy;
    }

    while (is.good())
    {
        std::string gene;
        std::string species;
        if (is >> gene)
        {
            if (is >> species)
            {
                gene2species.insert(gene, species);
            }
            else
            {
                std::ostringstream line_str;
                line_str << "Line " << lineno;
                is.close();
                throw AnError("The gene-to-species mapping seems to be "
                              "badly formatted. ", line_str.str());
            }
        }
        lineno++;
    }
    is.close();
    return gene2species;
}

// auwn to get the order on leaves
//////////////////////////////////////////////////////////
std::map<std::string, int>
TreeIO::readGeneSpeciesOrder(const std::string &filename)
{

    //std::vector<std::string> insertOrder; 
    std::map<std::string, int> io_map;
    
    int counter;	

    std::string last_species;

    ifstream is(filename.c_str());
    int lineno = 1;

    if(is.peek() == '#') // gs may start with a '#'
    {
        char dummy[LINELENGTH];
        is.getline(dummy, LINELENGTH);// >> dummy;
    }

    while (is.good())
    {
        std::string gene;
        std::string species;
        if (is >> gene)
        {
            if (is >> species)
            {
                if ( last_species.compare(species) != 0) 
                {
                  counter = 0;
                }
                
                last_species = species;
                io_map.insert(pair<string,int>(gene, ++counter));
                //insertOrder.push_back(gene);	   
            }
            else
            {
                std::ostringstream line_str;
                line_str << "Line " << lineno;
                is.close();
                throw AnError("The gene-to-species mapping seems to be "
                              "badly formatted. ", line_str.str());
            }
        }
        lineno++;
    }

     /*
     map<string, int>::iterator it;
     for ( it = io_map.begin(); it != io_map.end(); it++ )
     {
             std::cout << it->first  // string (key)
              << ':'
              << it->second          // string's value 
              << std::endl ;
      }
      */

    is.close();
    return io_map;
}
/////////////////////////////////////////////////////////

std::vector<StrStrMap>
TreeIO::readGeneSpeciesInfoVector(const std::string &filename)
{
    ifstream is(filename.c_str());
    int lineno = 1;
    std::vector<StrStrMap> gene2speciesVec;
    StrStrMap gene2species;
    std::string gsmark;
    is >> gsmark;
    if(gsmark != "#")
    {
        std::cerr << "error in gs vector, every gs must be preceeded by '#' line\n";
        exit(1);
    }

    while (is.good())
    {
        char dummy[LINELENGTH];
        is.getline(dummy, LINELENGTH);
        std::string gene;
        std::string species;
        if (is >> gene)
        {
            if(gene == "#")
            {
                gene2speciesVec.push_back(gene2species);
                gene2species.clearMap();
            }
            else
            {
                if (is >> species)
                {
                    gene2species.insert(gene, species);
                }
                else
                {
                    std::ostringstream line_str;
                    line_str << "(Line " << lineno << ")";
                    throw AnError("The gene-to-species mapping seems to be "
                                  "badly formatted. ", line_str.str());
                }
            }
        }
        
        lineno++;
    }
    gene2speciesVec.push_back(gene2species);
    return gene2speciesVec;
}

struct NHXtree *
        TreeIO::checkTagsForTree(TreeIOTraits& traits)
{
    // Reset all argmuents before starting
    traits.setNW(true);
    traits.setET(true);
    traits.setNT(true);
    traits.setBL(true);
    traits.setGS(true);
    traits.setAC(false);
    traits.setHY(false);

    // Parse tree from its source, into temporary datastructure
    struct NHXtree* t = readTree();
    struct NHXtree* ct = t;
    if (ct == 0)
    {
        throw AnError("The input tree is NULL!",
                      "TreeIO::checkTagsForTree()",
                      1);
    }

    while(ct)
    {
        if(TreeIO::recursivelyCheckTags(ct->root, traits) == false)
        {
            throw AnError("The input tree was empty!",
                          "TreeIO::checkTagsForTree()",
                          1);
        }
        ct= ct->next;
    }
    return t;
}


// Find the right value for edge time
Real
TreeIO::decideEdgeTime(struct NHXnode *v, const TreeIOTraits& traits,
                       bool isHY)
{
    Real edge_time = 0.0;

    if(traits.hasET()) // Use edge time info from file
    {
        if(traits.hasNWisET())
        {
            if(struct NHXannotation* a = find_annotation(v, "NW"))
            {
                edge_time = a->arg.t;
            }
            else if (isRoot(v))
            {
                edge_time = 0.0;
            }
            else
            {
                throw AnError("Edge without edge time found in tree.", 1);
            }
        }
        else if(struct NHXannotation *a = find_annotation(v, "ET"))
        {
            edge_time = a->arg.t;
        }
        else if (isRoot(v))
        {
            edge_time = 0.0;
        }
        else
        {
            throw AnError("Edge without edge time found in tree.", 1);
        }
        // Check for sanity
        if(edge_time <= 0)
        {
            if(edge_time < 0)
            {
                throw AnError("Tree contains an edge with negative time",1);
            }
            else if(isHY == false && !isRoot(v))
            {
                //TOFIX this is strange, this error is present only when using MAC
                throw AnError("Tree contains an edge with zero time.", 1);
            }
        }
    }

    return edge_time;
}

// Always include name , if it exists
std::string
TreeIO::decideNodeName(struct NHXnode *v) 
{
    std::string name = "";		// Default name is empty
    if (v->name)		// Leaf names are always inserted into 'name'
    {
        name = v->name;
    }
    else			// Inner nodes might have a name. Go check
    {
        struct NHXannotation *a = find_annotation(v, "S");
        if (a)			// Great, we found it
        {
            name = a->arg.str;	// Pick string from union 'arg'
        }
    }
    return name;
}

void
TreeIO::handleBranchLengths(Node *node, struct NHXnode *v, bool NWIsET)
{
    if(struct NHXannotation *a = find_annotation(v, "BL"))
    {
        node->setLength(a->arg.t);
    }
    else if(NWIsET)
    {
        throw AnError("TreeIO::extendBeepTree(...):\n"
                      "No branch length info found either in 'BL' and 'NW' is used for 'ET'",
                      234);
    }
    else if(struct NHXannotation *a = find_annotation(v, "NW"))
    {           // use info in 'NW'
        node->setLength(a->arg.t);
    }
    else if (v->parent)		// If not root...
    {
        throw AnError("TreeIO::extendBeepTree(...):\n"
                      "No branch length info found either in 'BL' or 'NW'",
                      234);
    }
}

// Recursively checks what tags are given for all nodes in subtree T_v
// Precondition: All bool argument has proper values. Assume a specific
// bool argument, 'A' has incoming value 'a', and the value for the 
// current subtree is 'b', then on return, A = a && b, i.e., false if 
// either a or b is false.
// postcondition: return true if subtree is non-empty, i.e, v != 0 
//----------------------------------------------------------------------
bool 
TreeIO::recursivelyCheckTags(struct NHXnode* v, TreeIOTraits& traits)
{
    if (v == 0) // i.e. if the parent was not a leaf
    {
        return false;
    }

    recursivelyCheckTags(v->left, traits);
    recursivelyCheckTags(v->right, traits);
    checkTags(*v, traits);
    return true;
}


// Checks what tags are given for node v
// Precondition: All bool argument has proper values. Assume a specific
// bool argument, 'A' has incoming value 'a', and the value for the 
// current node is 'b', then on return, A = a && b.
//----------------------------------------------------------------------
void
TreeIO::checkTags(struct NHXnode& v, TreeIOTraits& traits)
{
    // Determine if NW is given
    //-------------------------------------------------------------------
    if(!find_annotation(&v, "NW") && !isRoot(&v))
    {
        traits.setNW(false);
    }

    // Determine if ET is given
    //-------------------------------------------------------------------
    if(!find_annotation(&v, "ET") && !isRoot(&v))
    {
        traits.setET(false);
    }

    // Check if NT is given
    //-------------------------------------------------------------------
    if(!find_annotation(&v, "NT") && !isLeaf(&v))
    {
        traits.setNT(false);
    }

    // Check if BL is given
    //-------------------------------------------------------------------
    if(!find_annotation(&v, "BL") && !isRoot(&v))
    {
        traits.setBL(false);
    }

    // Check if AC is given.
    //! \todo{The AC check behaves conversely from other checks. Any presence of
    //! AC makes the hasAC argument true.}
    //-------------------------------------------------------------------
    if(find_annotation(&v, "AC"))
    {
        traits.setAC(true);
    }

    // Check if GS is given for leaves.
    //-------------------------------------------------------------------
    if (v.left == 0 && v.right == 0 && speciesName(&v) == 0)
    {
        traits.setGS(false);
    }

    // Check if there are hybrid annotations
    if(find_annotation(&v, "HY") || find_annotation(&v, "EX") ||
            find_annotation(&v, "OP"))
    {
        traits.setHY(true);
    }
}


// Generic reading function, interfacing NHX* code
//----------------------------------------------------------------------
NHXtree*
TreeIO::readTree()
{
    if (source == readFromStdin)
    {
        return read_tree(NULL);
    }
    else if (source == readFromFile)
    {
        return read_tree(stringThatWasPreviouslyNamedS.c_str());
    }
    else if (source == readFromString)
    {
        return read_tree_string(stringThatWasPreviouslyNamedS.c_str());
    }
    else
    {
        throw AnError("TreeIO not properly initialized!");
        return NULL;
    }
}
