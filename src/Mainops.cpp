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

#include "Mainops.h"

#include "tree/Node.hh"
#include "Parameters.h"
#include "tree/TreeIO.hh"
#include "tree/Node.hh"
#include "tree/TreeIO.hh"
#include "utils/AnError.hh"
#include "draw/DrawTreeCairo.h"
#include "layout/Layoutrees.h"

#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp> 

using namespace std;

Mainops::Mainops()
    :Guest(0),Host(0),gamma(0),lambdamap(0),dt(0),parameters(0),io(0)
{
	std::cerr << "1) Mainops constructor called" << std::endl;
}

void Mainops::start()
{
    std::cerr << "2) Mainops start called" << std::endl;
    io = new TreeIO();
    dt = new DrawTreeCairo();
}

void Mainops::cleanTrees()
{
    if(Guest)
    {
        delete(Guest);
    }
    Guest = 0;
    

    if(Host)
    {
        delete(Host);
    }
    Host = 0;

    if(io)
    {
        delete(io);
    }
    io = 0;

    io = new TreeIO();
    AC.clear();
    gs.clearMap();
}

Mainops::~Mainops()
{
   
    if(Guest)
    {
        delete(Guest);
    }
    Guest = 0;

    if(Host)
    {
        delete(Host);
    }
    Host = 0;

    if(gamma)
    {
        delete(gamma);

    }
    gamma = 0;

    //if(lambdamap)
    //{
        //delete(lambdamap);

    //}
    //lambdamap = 0;

    if(dt)
    {
        delete(dt);

    }
    dt = 0;

    if(io)
    {
        delete(io);
    }
    io = 0;
}

bool Mainops::lateralTransfer(const std::string &mapname, bool dp)
{
    std::cerr << "XXX Mainops loadPreComputedScenariolateralTransfer called" << std::endl;	

    Phyltr late = Phyltr();
    late.g_input.duplication_cost = parameters->lateralduplicost;
    late.g_input.transfer_cost = parameters->lateraltrancost;
    late.g_input.max_cost = parameters->lateralmaxcost;
    late.g_input.min_cost = parameters->lateralmincost;
    late.g_input.gene_tree = Guest;
    late.g_input.species_tree = Host;
    if(dp)
    {
        late.g_input.print_only_minimal_loss_scenarios = false;
        late.g_input.print_only_minimal_transfer_scenarios = false;
    }
    if(parameters->isreconciled)
    {
        late.g_input.sigma_fname = mapname;
        late.read_sigma();
    }
    else
    {
        late.read_sigma(gs.getMapping());
    }
    
    if(dp)
    {
        late.dp_algorithm();
        late.backtrack();
    }
    else
    {
        late.fpt_algorithm();
    }

    if(late.scenarios.size() > 0 && thereAreLGT(late.scenarios))
    {
        Scenario scenario = late.getMinCostScenario();
        lambda = scenario.cp.getLambda();
        transferedges = scenario.transfer_edges;
        parameters->transferedges = transferedges;
        sigma = late.g_input.sigma;
        scenarios = late.scenarios;
        return true;
    }
    else
    {
        parameters->lattransfer = false;
        std::cerr << "Not valid LGT scenarios found.." << std::endl;
        return false;
    }
}

void Mainops::printLGT()
{
    std::cerr << "List of computed LGT scenarios sorted by cost.." << std::endl;
    BOOST_FOREACH(Scenario &sc, scenarios)
    {
        std::cerr << sc << std::endl;
    }
}

bool Mainops::thereAreLGT(const std::vector<Scenario> &scenarios) const
{
     std::cerr << "XXX Mainops thereAreLGT called" << std::endl;	
   BOOST_FOREACH(const Scenario &sc, scenarios)
    {
        if(sc.transfer_edges.any())
        {
            return true;
        }
    }
    return false;
}

void Mainops::OpenReconciled(const string &reconciled)
{
     std::cerr << "XXX Mainops OpenReconciled called" << std::endl;	
    io->setSourceFile(reconciled);
    Guest = new TreeExtended(io->readBeepTree<TreeExtended,Node>(&AC, &gs));
}

void Mainops::OpenHost(const string &species)
{    
     std::cerr << "XXX Mainops OpenHost called" << std::endl;	
    io->setSourceFile(species);
    Host = new TreeExtended(io->readHostTree<TreeExtended,Node>());
    Node *root = Host->getRootNode();

    if ((double)root->getTime() != (double)0.0)
    {
        Real t = root->getNodeTime();
        root->setTime(0.1 * t);
    }

    if (Host->imbalance() / Host->rootToLeafTime() > 0.01)
    {
        parameters->scaleByTime = false;
        std::cerr << "The species tree is not ultrametric (it appears unbalanced),\n"
                     "so scaling by time is turned off. See also option '-t'" << std::endl;

    }
}

void Mainops::CalculateGamma()
{
     std::cerr << "XXX CalculateGamma OpenHost called" << std::endl;	
  
    if (parameters->do_not_draw_species_tree == false)
    {
        gamma = new GammaMapEx<Node>(*Guest, *Host, gs, AC);
        if (parameters->lattransfer)
        {
            gamma = new GammaMapEx<Node>(gamma->update(*Guest,*Host,sigma,transferedges));
        }
        lambdamap = gamma->getLambda();
    }
    else
    {
        lambdamap = new LambdaMapEx<Node>(*Guest, *Host, gs);
        if (parameters->lattransfer)
        {
            lambdamap->update(*Guest,*Host,sigma,transferedges);
        }
        gamma = new GammaMapEx<Node>(GammaMapEx<Node>::MostParsimonious(*Guest,*Host,*lambdamap));
    }
}

void Mainops::reconcileTrees(const string &gene, const string &species, const string &mapfile)
{
    std::cerr << "3) Mainops reconcileTrees called" << std::endl;	
    io->setSourceFile(gene);
    Guest = new TreeExtended(io->readBeepTree<TreeExtended,Node>(&AC, &gs));
    io->setSourceFile(species);
    Host = new TreeExtended(io->readNewickTree<TreeExtended,Node>());

    if(mapfile != "")
    {
        gs = TreeIO::readGeneSpeciesInfo(mapfile);
    }
    else
    {
        throw AnError(": The mapfile is empty!\n");
    }
    
    
    /*
    for(std::vector<std::string>::iterator it = insertOrder.begin(); it != insertOrder.end(); ++it) {
       std::cerr << *it << endl;
    }
    */
    

    // Print the ordered map
    /*
    std::vector<std::string> insertOrder =  gs.getInsertionOrder();
    std::cerr << insertOrder.size()<< '\n';
    for (int i = 0; i < insertOrder.size(); ++i)
    {
       const std::string &s = insertOrder[i];
       std::cerr << s << ' ' << gs.find(s) << '\n';
    }	
    */

   /*
     std::map<std::string, std::string> gsmap =  getMapping();
    for (map<string,string>::const_iterator i = gsmap.begin(); i != gsmap.end(); i++)
    {

             std::cerr << i->first + "\t" + i->second + "\n" << std::endl;

    }
    */
    

    LambdaMapEx<Node> lambdamap_local = LambdaMapEx<Node>(*Guest, *Host, gs);
    GammaMapEx<Node> gamma_local = GammaMapEx<Node>(GammaMapEx<Node>::MostParsimonious(*Guest, *Host, lambdamap_local));
    string textTree = io->writeGuestTree<TreeExtended,Node>(*Guest,&gamma_local);
    io->setSourceString(textTree);
    Guest = new TreeExtended(io->readBeepTree<TreeExtended,Node>(&AC, &gs));


    std::map<std::string, int> io_map  =  TreeIO::readGeneSpeciesOrder(mapfile);
    /*
    // To show the leaf order
    std::cerr <<"Hello size = " << io_map.size() << std::endl; 
    map<string, int>::iterator it;
    for ( it = io_map.begin(); it != io_map.end(); it++ )
    {
          std::cerr << it->first  // string (key)
          << ':'
          << it->second          // string's value 
          << std::endl ;
     }
    */


    int position;
    map<string, int>::iterator iter;
    for (unsigned i = 0; i < Guest->getNumberOfNodes(); i++) {
     Node *n   = Guest->getNode(i);
     if(n->isLeaf()){
	const string &leafName = n->getName();
        iter = io_map.find(leafName);
        if (iter != io_map.end())
        position = iter->second;
     	n->setPosition(position-1);		
     }	    
    }
    io_map.clear(); 
    


    OpenHost(species);
}

void Mainops::calculateCordinates()
{
    std::cerr << "XXX Mainops calculateCordinates called" << std::endl;	
    Host->reset();
    Guest->reset();
    //reduce crossing only if not LGT
    if(parameters->reduce && !(bool)(parameters->lattransfer))
    {
        std::cerr << "NOTE : the option -R is still experimental.." << std::endl;
        gamma->twistAndTurn();
    }
    LayoutTrees spcord = LayoutTrees(Host,Guest,parameters,gamma,lambdamap);
    spcord.start();
    parameters->leafwidth = spcord.getNodeHeight();
    //printf("%f",parameters->leafwidth);
}

int Mainops::checkValidity()
{
    std::cerr << "XXX Mainops checkValidity called" << std::endl;	
    if(parameters->lattransfer)
    {
        return getValidityLGT();
    }
    else
    {
        return gamma->valid();
    }
}

void Mainops::DrawTree(cairo_t *cr)
{
    std::cerr << "5) Mainops DrawTree called" << std::endl;	

    bool heatMapMode = true;
    dt->start(parameters,Guest,Host,gamma,lambdamap,cr);
    dt->setHeatMap(heatMapMode);	
    dt->calculateTransformation();

    
    if(parameters->do_not_draw_species_tree == false)
    {
        //species tree
        if (!parameters->noTimeAnnotation)
        {
            if (parameters->timeAtEdges)
            {
                dt->TimeLabelsOnEdges();
            }
            else
            {
                dt->DrawTimeEdges();
                dt->DrawTimeLabels();
            }
        }
        //dt->DrawSpeciesEdgesWithContour(); //NOTE this breaks the layout
        dt->DrawSpeciesEdges();
        dt->DrawSpeciesNodes();
        dt->DrawSpeciesNodeLabels();
    }
    // gene tree
    if(!parameters->do_not_draw_guest_tree)
    {
        dt->DrawGeneEdges();
        dt->DrawGeneNodes();
        dt->DrawGeneLabels();  //auwn: dont draw gene labels
        if(parameters->markers)
        {
            dt->GeneTreeMarkers();
        }
    }
    
    if(parameters->header)
    {
        dt->createHeader();
    }
    if(parameters->legend)
    {
        dt->createLegend();
    }
    if(parameters->title)
    {
        dt->createTitle();
    }
    if(parameters->show_event_count)
    {
        dt->writeEventCosts();
    }


    // Legend for heatMap here ...
    if(heatMapMode)
      dt->createMyLegend();
    
}

int Mainops::RenderImage()
{
    std::cerr << "6) Mainops RenderImage called" << std::endl;	
  
    bool ok = dt->RenderImage();
    dt->cleanUp();
    return ok;
}

Parameters* Mainops::getParameters()
{
  
    return parameters;
}

void Mainops::setParameters(Parameters *p)
{
    parameters = p;
}

bool Mainops::getValidityLGT()
{
    if(gamma->validLGT())
    {
        return true;
    }
    else
    {
        sort(scenarios.begin(), scenarios.end());

        BOOST_FOREACH (Scenario &sc, scenarios)
        {
            transferedges = sc.transfer_edges;
            parameters->transferedges = sc.transfer_edges;
            parameters->duplications = sc.duplications;
            lambda = sc.cp.getLambda();
            CalculateGamma();
            if(gamma->validLGT())
            {
                return true;
            }
        }
        return false;
    }
}

void Mainops::drawBest()
{
    std::cerr << "4) Mainops drawBest called" << std::endl;	
    
    CalculateGamma(); //calculation of gamma and lambda
    if(!checkValidity())
    {
        if(parameters->lattransfer)
        {
            throw AnError(": The LGT scenario was not valid. Aborts!\n");
        }
        else
        {
            throw AnError(": This is not a correctly reconciled tree. Aborts!\n");
        }
    }
    calculateCordinates(); //calculation of the drawing cordinates
    DrawTree();  //drawing the tree
    if(!RenderImage())
      {
      	throw AnError(": Error rendering the tree..\n");
      } // save the file
}

void Mainops::drawAllLGT()
{
      std::cerr << "XXX Mainops drawAllLGT called" << std::endl;	

    unsigned index = 0;
    std::string original_filename = parameters->outfile;
    sort(scenarios.begin(), scenarios.end());
    BOOST_FOREACH (Scenario &sc, scenarios)
    {
        transferedges = sc.transfer_edges;
        parameters->transferedges = sc.transfer_edges;
        parameters->duplications = sc.duplications;
        lambda = sc.cp.getLambda();
        parameters->outfile = original_filename + boost::lexical_cast<string>(++index);
        CalculateGamma(); //calculation of gamma and lambdamap
        if(gamma->validLGT())
        {
            calculateCordinates(); //calculation of the drawing cordinates
            DrawTree();  //drawing the tree
            if(!RenderImage())
            {
                throw AnError(": Error rendering the tree..\n");
            } // save the file
        }
    }
}

void Mainops::loadPreComputedScenario(const std::string &filename,const std::string &mapname)
{
  std::cerr << "XXX Mainops loadPreComputedScenario called" << std::endl;	

    std::ifstream scenario_file;
    std::string line;
    scenario_file.open(filename.c_str(), std::ios::in);

    if (!scenario_file)
    {
        throw AnError("Could not open file " + filename);
    }

    while (getline(scenario_file, line))
    {
        if(scenario_file.good())
        {
            if(line.size() == 0)  // Skip any blank lines
            {
                continue;
            }
            else if(line[0] == '#')  // Skip any comment lines
            {
                continue;
            }
            else if ((line.find("Transfer") != std::string::npos))
            {
                const std::size_t start_pos = line.find(":");
                const std::size_t stop_pos = line.size() - 1;
                std::string temp = line.substr(start_pos + 1,stop_pos - start_pos);
                stringstream lineStream(temp);
                std::vector<std::string> transfer_nodes((istream_iterator<std::string>(lineStream)), istream_iterator<std::string>());
                transferedges.clear();
                transferedges.resize(Guest->getNumberOfNodes());
                for(std::vector<std::string>::const_iterator it = transfer_nodes.begin(); it != transfer_nodes.end(); ++it)
                {
                    std::vector<std::string> strs;
                    std::string temp = *it;
                    temp.erase(remove(temp.begin(),temp.end(),'('),temp.end());
                    temp.erase(remove(temp.begin(),temp.end(),')'),temp.end());
                    boost::split(strs, temp, boost::is_any_of(","));
                    if(Guest->getNode(boost::lexical_cast<unsigned>(strs.at(0))) != NULL)
                    {
                        transferedges.set(boost::lexical_cast<unsigned>(strs.at(0))); //origin LGT
                    }
                    else
                    {
                        throw AnError("Node read in the LGT scenario file does not exist in the Gene Tree");
                    }
                    //transferedges.set(boost::lexical_cast<unsigned>(strs.at(1))); //destiny LGT
                    // strs.at(2) //this is the time NOTE not used yet (the idea is to put the time in the node and use it to compute the cordinates
                }
            }
        }
    }
    scenario_file.close();
    parameters->lattransfer = true;
    Phyltr late = Phyltr();
    late.g_input.gene_tree = Guest;
    late.g_input.species_tree = Host;

    if(parameters->isreconciled)
    {
        late.g_input.sigma_fname = mapname;
        late.read_sigma();
    }
    else
    {
        late.read_sigma(gs.getMapping());
    }

    parameters->transferedges = transferedges;
    sigma = late.g_input.sigma;
    return;
}
