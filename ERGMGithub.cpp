//============================================================================
// Name        : SEPNET
// Author      : Hazem Krichene
// Copyright   : PostK project - University of Hyogo
// Description : Simulation and Estimation of Production Network
//============================================================================

#include <iostream>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <iterator>
#include <set>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <numeric>
#include <math.h>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
//==============================================================================

struct Link{
	std::set<std::string> in_neighbors; /*Contains the list of suppliers of firm ID*/

	std::set<std::string> out_neighbors; /*Contains the list of customers of firm ID*/
};

//==============================================================================
int number_variables = 16;
int IFD_size = 2; // The size of the MCMC sampling
VectorXd theta(number_variables);
int initial_simulation = 5; //Initialization to calculate Theta
double a_r = 0.1; //Initial "a" parameter: See Lusher et al. (2013)
int number_phases = 10; //Define the maximum phases until convergence; Value of 10 means that if convergence does not occur after 10 phases, the code is considered as non converging
int burnin = 100; //The burnin MCMC is used to avoid the effect to initialization
int NumNetworks = 3; //How many network do you want to simulate after the convergence to achieve the GoF analysis
//==============================================================================

/*
 * network: key is firm ID, with sets of suppliers and customers
 * attribute: key is firm ID, with vector of attributes (location, industrial sector, etc.)
 * GPS_Distance: Log-Distance between prefecture P1 and prefecture P2
 * Firms_Sectors: Supplier ID with the list of industrial sectors of its customers
*/

std::unordered_map<std::string, struct Link> network;
std::unordered_map<std::string, vector< std::string> > attribute;
unordered_map< std::string, unordered_map< std::string, double> > GPS_Distance;
unordered_map< std::string, set<std::string>> Firms_Sectors;

//===============================Function to Upload data ============================================================================

void Upload_Data()
{
	ifstream attribute_file, prod, gps;

	//===========================================================================================

	/*
	 * Files structures to compile the program:
	 * Text file Attributes contains 6 Columns: Firm ID, Prefecture code, Industry Code, Major Bank Code, Number of Employees, Total Sales
	 * Text file Edges contains 2 columns: Supplier ID, Customer ID
	 * Text file Geo_Distance contains 3 columns: Prefecture I Code, Prefecture J Code, Log-Distance between prefectures I and J
	*/

	attribute_file.open("ToyAttributes.dat");
	prod.open("ToyEdges.dat");
	gps.open("Geo_Distance(Real Japanese Prefectures).dat");


	//============================================================================================

	std::string line_att, line_prod, line_gps;

	while (getline(attribute_file,line_att))
	{
	    istringstream ss(line_att);

	    std::string FirmID;
	    std::string PrefCode;
	    std::string IndustryCode;
	    std::string BankID;
	    std::string Employee;
	    std::string Sales;

		ss >> FirmID;
		ss >> PrefCode;
		ss >> IndustryCode;
		ss >> BankID;
		ss >> Employee;
		ss >> Sales;

		attribute[FirmID].push_back(PrefCode);
        attribute[FirmID].push_back(IndustryCode);
        attribute[FirmID].push_back(BankID);
        attribute[FirmID].push_back(Employee);
        attribute[FirmID].push_back(Sales);
	}

	while(getline(prod,line_prod))
	{
		std::string supplier;
		std::string customer;

		istringstream ss(line_prod);
        ss >> supplier;
	    ss >> customer;
        network[supplier].out_neighbors.insert(customer);
    	network[customer].in_neighbors.insert(supplier);

    	Firms_Sectors[supplier].insert(attribute[customer][1]);
	}

	prod.close();

	while(getline(gps,line_gps))
	{
		std::string pref1;
		std::string pref2;
		std::string distance;

		istringstream ss(line_gps);
        ss >> pref1;
	    ss >> pref2;
	    ss >> distance;

    	GPS_Distance[pref1][pref2] = std::stod(distance);
	}

	gps.close();
}
 // ========================================= Network statistics functions =================================================================================

std::vector<std::string> Set_Intersection(std::set<std::string>& set1, std::set<std::string>& set2)
{
	std::vector<std::string> v_intersection;

    std::set_intersection(set1.begin(), set1.end(),
    		set2.begin(), set2.end(),std::back_inserter(v_intersection));

    return v_intersection;
}

std::set<std::string> Neighbors(std::string FirmID, std::string& Mode)
{
	if (Mode == "IN")
	{
		return network[FirmID].in_neighbors;
	}

	else if (Mode == "OUT")
	{
		return network[FirmID].out_neighbors;
	}

	else if (Mode == "ALL")
	{
		std::set<std::string> all_neighbors = network[FirmID].in_neighbors;
		all_neighbors.insert(network[FirmID].out_neighbors.begin(), network[FirmID].out_neighbors.end());
		return all_neighbors;
	}
}

bool Adjacency(std::string Firm_i, std::string Firm_j)
{
	return network[Firm_i].out_neighbors.find(Firm_j) != network[Firm_i].out_neighbors.end();
}

vector<int> Degree()
{
    int d = 0;
    int r = 0;
    for (std::unordered_map<std::string, struct Link>::iterator itr=network.begin(); itr!=network.end();itr++)
    {
    	std::set<std::string> out_neighbor  = network[(*itr).first].out_neighbors;
        d += out_neighbor.size();
        for (std::set<std::string>::iterator itr_out=out_neighbor.begin(); itr_out!=out_neighbor.end();++itr_out)
        {
        	if (Adjacency(*itr_out,(*itr).first) == 1) r++ ;
        }
    }
    vector <int> v;
    v.push_back(d);
    v.push_back(r/2);
    return v;
}

int InDegree(std::string& FirmID)
{
	return network[FirmID].in_neighbors.size();
}

int OutDegree(std::string& FirmID)
{
	return network[FirmID].out_neighbors.size();
}

vector <int> InDegree_Distribution()
{
    vector <int> in_degree_distribution;
    for (std::unordered_map<std::string, struct Link>::iterator itr=network.begin(); itr!=network.end();itr++)
    {
    	in_degree_distribution.push_back(network[(*itr).first].in_neighbors.size());
    }
	return in_degree_distribution;
}

vector <int> OutDegree_Distribution()
{
    vector <int> out_degree_distribution;
    for (std::unordered_map<std::string, struct Link>::iterator itr=network.begin(); itr!=network.end();itr++)
    {
    	out_degree_distribution.push_back(network[(*itr).first].out_neighbors.size());
    }
	return out_degree_distribution;
}

// ============================================== End of Network Statistics functions =================================================================

// ============================================== Calculate the global network configuration =================================================================

vector<double> Network_Configuration()
{
	int z_edges = Degree()[0];
    int z_reciprocity = Degree()[1];
    std::string Mode_Out = "OUT";
    std::string Mode_In = "IN";

    vector<double> z_triangles;
    z_triangles.push_back(0.0);
    z_triangles.push_back(0.0);
    z_triangles.push_back(0.0);
    z_triangles.push_back(0.0);
    double z_two_path = 0.0;
    double profit_heterophily = 0.0;
    double profit_sender = 0.0;
    double profit_receiver = 0.0;
    int sector_homophily = z_edges;
    int bank_homophily = 0;
    double geographic_homophily = 0.0;
    int employee_heterophily = 0;
    for (std::unordered_map<std::string, struct Link>::iterator source_node=network.begin(); source_node!=network.end();++source_node)
    {
    	for (std::set<std::string>::iterator target_node=network[(*source_node).first].out_neighbors.begin(); target_node!=network[(*source_node).first].out_neighbors.end();++target_node)
    	{
    		std::set<std::string> out_neighbor_source = Neighbors((*source_node).first,Mode_Out);
		    std::set<std::string> in_neighbor_source = Neighbors((*source_node).first,Mode_In);
		    std::set<std::string> out_neighbor_target = Neighbors(*target_node,Mode_Out);
		    std::set<std::string> in_neighbor_target = Neighbors(*target_node,Mode_In);

		    int in_out = Set_Intersection(in_neighbor_source,out_neighbor_target).size();
		    int out_in = Set_Intersection(out_neighbor_source,in_neighbor_target).size();
		    int out_out = Set_Intersection(out_neighbor_source,out_neighbor_target).size();
		    int in_in = Set_Intersection(in_neighbor_source,in_neighbor_target).size();

		    z_triangles[0] += Adjacency((*source_node).first,*target_node)*2*(1 - pow(0.5,in_out));
		    z_triangles[1] += Adjacency((*source_node).first,*target_node)*2*(1 - pow(0.5,out_in));
		    z_triangles[2] += Adjacency((*source_node).first,*target_node)*2*(1 - pow(0.5,out_out));
		    z_triangles[3] += Adjacency((*source_node).first,*target_node)*2*(1 - pow(0.5,in_in));

		    z_two_path += 2*(1 - pow(0.5,in_out)) + 2*(1 - pow(0.5,out_in));

		    employee_heterophily += abs(atof(attribute[(*source_node).first][3].c_str()) - atof(attribute[*target_node][3].c_str()));

		    profit_heterophily += abs(atof(attribute[(*source_node).first][4].c_str()) - atof(attribute[*target_node][4].c_str()));
		    profit_sender += atof(attribute[(*source_node).first][4].c_str());
		    profit_receiver += atof(attribute[*target_node][4].c_str());
		    if (attribute[(*source_node).first][2] == attribute[*target_node][2]) bank_homophily++;
		    geographic_homophily+=GPS_Distance[attribute[(*source_node).first][0]][attribute[*target_node][0]];

    	}
    }

    vector <int> out_degrees = OutDegree_Distribution();
    int max_out = *std::max_element( out_degrees.begin(), out_degrees.end());
    std::vector<int> out_dist(max_out+1, 0);
    for (std::vector<int>::iterator  itr = out_degrees.begin(); itr!= out_degrees.end(); ++itr) out_dist[*itr]++;
    std::vector<double> GWOD;
    for (std::vector<int>::iterator  itr = out_dist.begin(); itr!= out_dist.end(); ++itr)
    {
    	int position = std::distance(out_dist.begin(), itr);
    	GWOD.push_back(exp(-(double)position)*(double)(*itr));
    }
    double z_out_stars = std::accumulate(GWOD.begin(), GWOD.end(), 0);

    vector <int> in_degrees = InDegree_Distribution();
    int max_in = *std::max_element( in_degrees.begin(), in_degrees.end());
    std::vector<int> in_dist(max_in+1, 0);
    for (std::vector<int>::iterator  itr = in_degrees.begin(); itr!= in_degrees.end(); ++itr) in_dist[*itr]++;
    std::vector<double> GWID;
    for (std::vector<int>::iterator  itr = in_dist.begin(); itr!= in_dist.end(); ++itr)
    {
    	int position = std::distance(in_dist.begin(), itr);
    	GWID.push_back(exp(-(double)position)*(double)(*itr));
    }
    double z_in_stars = std::accumulate(GWID.begin(), GWID.end(), 0);
    vector<double> configuration;
    configuration.push_back(z_edges);
    configuration.push_back(z_reciprocity);
    configuration.push_back(sector_homophily);
    configuration.push_back(-geographic_homophily);
    configuration.push_back(profit_heterophily);
    configuration.push_back(profit_sender);
    configuration.push_back(profit_receiver);
    configuration.push_back(z_in_stars);
    configuration.push_back(z_out_stars);
    configuration.push_back(z_two_path);
    configuration.push_back(z_triangles[0]);
    configuration.push_back(z_triangles[1]);
    configuration.push_back(z_triangles[2]);
    configuration.push_back(z_triangles[3]);
    configuration.push_back(bank_homophily);
    configuration.push_back(employee_heterophily);

    return configuration;
}

// ============================================== End of Calculation of the global network configuration =================================================================

// ============================================== Functions to calculate the change statistics ============================================================================

/*
 * Functions related to endogenous network statistics like K_Triangles, K_Stars, etc., may be used for any work.
 * Functions related to exogenous network statistics are adapted to the production network with the special case of our attributes.
 * See Krichene et al. (2018) for further details.
*/

int Edge_Change(int xij)
{
	if (xij == 1)
	{
		return 1;
	}

	else
	{
		return -1;
	}
}


int Reciprocity_Change(std::string node_i, std::string node_j,int xij)
{
    int z_reciprocity = 0;
    if ((xij == 1) && Adjacency(node_j, node_i) == 1)
	{
    	z_reciprocity = 1;
	}

    else if ((xij == 0) && Adjacency(node_j, node_i) == 1)
	{
    	z_reciprocity = -1;
	}

    else z_reciprocity = 0;
    return z_reciprocity;
}


int Sector_Change(std::string node_i, std::string node_j,int xij)
{
    int sector_indicator = 0;
    int z_sector = 0;

    if (Firms_Sectors[node_i].find(attribute[node_j][1]) != Firms_Sectors[node_i].end()==1)
    	sector_indicator = 1;

    if (xij == 1)
        z_sector = sector_indicator;
    else if (xij == 0)
        z_sector = -sector_indicator;

    return z_sector;
}

int Bank_Change(std::string node_i, std::string node_j,int xij)
{
    int bank_indicator = 0;
    int z_bank = 0;

    if (attribute[node_i][2] == attribute[node_j][2])
    	bank_indicator = 1;

    if (xij == 1)
    	z_bank = bank_indicator;
    else if (xij == 0)
    	z_bank = -bank_indicator;

    return z_bank;
}


double Geographic_Change(std::string node_i, std::string node_j,int xij)
{
    return -GPS_Distance[attribute[node_i][0]][attribute[node_j][0]];
}


vector<double> Profit_Change(std::string node_i, std::string node_j,int xij)
{
	vector<double> z_profit;
	z_profit.push_back(0.0);
	z_profit.push_back(0.0);
	z_profit.push_back(0.0);

    if (xij == 1)
    {
        z_profit[0] = abs(atof(attribute[node_i][4].c_str()) - atof(attribute[node_j][4].c_str()));
        z_profit[1] = atof(attribute[node_i][4].c_str());
        z_profit[2] = atof(attribute[node_j][4].c_str());
    }

    else if (xij == 0)
    {
        z_profit[0] = -abs(atof(attribute[node_i][4].c_str()) - atof(attribute[node_j][4].c_str()));
        z_profit[1] = -atof(attribute[node_i][4].c_str());
        z_profit[2] = -atof(attribute[node_j][4].c_str());
    }

    return z_profit;
}

int Employee_Change(std::string node_i, std::string node_j,int xij)
{
	int z_employee;

    if (xij == 1)
    {
    	z_employee = abs(atof(attribute[node_i][3].c_str()) - atof(attribute[node_j][3].c_str()));
    }

    else if (xij == 0)
    {
    	z_employee = -abs(atof(attribute[node_i][3].c_str()) - atof(attribute[node_j][3].c_str()));
    }

    return z_employee;
}

vector<double> Triangles_Change(std::string node_i, std::string node_j,int xij)
{

    if (xij == 1)
    {
		network[node_i].out_neighbors.erase(node_j);
		network[node_j].in_neighbors.erase(node_i);
    }
    else
    {
		network[node_i].out_neighbors.insert(node_j);
		network[node_j].in_neighbors.insert(node_i);
    }

    std::string Mode_Out = "OUT";
    std::string Mode_In = "IN";
    std::string Mode_All = "ALL";
    std::set<std::string> affected_nodes;
    affected_nodes.insert(node_i);
    affected_nodes.insert(node_j);
    std::set<std::string> all_neighbor_i = Neighbors(node_i,Mode_All);
    std::set<std::string> all_neighbor_j = Neighbors(node_j,Mode_All);
    std::vector<std::string> common_neighbors = Set_Intersection(all_neighbor_i,all_neighbor_j);
    affected_nodes.insert(common_neighbors.begin(), common_neighbors.end());

    vector<double> initial_triangle;
    initial_triangle.push_back(0.0);
    initial_triangle.push_back(0.0);
    initial_triangle.push_back(0.0);
    initial_triangle.push_back(0.0);
    vector<double> final_triangle;
    final_triangle.push_back(0.0);
    final_triangle.push_back(0.0);
    final_triangle.push_back(0.0);
    final_triangle.push_back(0.0);
    double initial_two_path = 0;
    double final_two_path = 0;

    for (std::set<std::string>::iterator source_node= affected_nodes.begin(); source_node!= affected_nodes.end();++source_node)
    {
    	for (std::set<std::string>::iterator target_node = network[*source_node].out_neighbors.begin();target_node != network[*source_node].out_neighbors.end();++target_node)
    	{
    		std::set<std::string> out_neighbor_source = Neighbors(*source_node,Mode_Out);
		    std::set<std::string> in_neighbor_source = Neighbors(*source_node,Mode_In);
		    std::set<std::string> out_neighbor_target = Neighbors(*target_node,Mode_Out);
		    std::set<std::string> in_neighbor_target = Neighbors(*target_node,Mode_In);

		    int in_out = Set_Intersection(in_neighbor_source,out_neighbor_target).size();
		    int out_in = Set_Intersection(out_neighbor_source,in_neighbor_target).size();
		    int out_out = Set_Intersection(out_neighbor_source,out_neighbor_target).size();
		    int in_in = Set_Intersection(in_neighbor_source,in_neighbor_target).size();

		    initial_triangle[0] += Adjacency(*source_node,*target_node)*2*(1 - pow(0.5,in_out));
		    initial_triangle[1] += Adjacency(*source_node,*target_node)*2*(1 - pow(0.5,out_in));
		    initial_triangle[2] += Adjacency(*source_node,*target_node)*2*(1 - pow(0.5,out_out));
		    initial_triangle[3] += Adjacency(*source_node,*target_node)*2*(1 - pow(0.5,in_in));

		    initial_two_path += 2*(1 - pow(0.5,in_out)) + 2*(1 - pow(0.5,out_in));
    	}
    }

    if (xij == 1)
    {
		network[node_i].out_neighbors.insert(node_j);
		network[node_j].in_neighbors.insert(node_i);
    }
    else
    {
		network[node_i].out_neighbors.erase(node_j);
		network[node_j].in_neighbors.erase(node_i);
    }
    for (std::set<std::string>::iterator source_node= affected_nodes.begin(); source_node!= affected_nodes.end();++source_node)
    {
    	for (std::set<std::string>::iterator target_node = network[*source_node].out_neighbors.begin();target_node != network[*source_node].out_neighbors.end();++target_node)
    	{
    		std::set<std::string> out_neighbor_source = Neighbors(*source_node,Mode_Out);
		    std::set<std::string> in_neighbor_source = Neighbors(*source_node,Mode_In);
		    std::set<std::string> out_neighbor_target = Neighbors(*target_node,Mode_Out);
		    std::set<std::string> in_neighbor_target = Neighbors(*target_node,Mode_In);

		    int in_out = Set_Intersection(in_neighbor_source,out_neighbor_target).size();
		    int out_in = Set_Intersection(out_neighbor_source,in_neighbor_target).size();
		    int out_out = Set_Intersection(out_neighbor_source,out_neighbor_target).size();
		    int in_in = Set_Intersection(in_neighbor_source,in_neighbor_target).size();

		    final_triangle[0] += Adjacency(*source_node,*target_node)*2*(1 - pow(0.5,in_out));
		    final_triangle[1] += Adjacency(*source_node,*target_node)*2*(1 - pow(0.5,out_in));
		    final_triangle[2] += Adjacency(*source_node,*target_node)*2*(1 - pow(0.5,out_out));
		    final_triangle[3] += Adjacency(*source_node,*target_node)*2*(1 - pow(0.5,in_in));

		    final_two_path += 2*(1 - pow(0.5,in_out)) + 2*(1 - pow(0.5,out_in));
    	}
    }

    vector<double> z_triangle;
    z_triangle.push_back(final_triangle[0] - initial_triangle[0]);
    z_triangle.push_back(final_triangle[1] - initial_triangle[1]);
    z_triangle.push_back(final_triangle[2] - initial_triangle[2]);
    z_triangle.push_back(final_triangle[3] - initial_triangle[3]);
    z_triangle.push_back(final_two_path - initial_two_path);
    return z_triangle;
}


vector<double> Stars_Change(std::string node_i, std::string node_j,int xij)
{
    int ind;
    if (xij == 0) ind = 1;
    else ind = -1;

    int d_in = InDegree(node_j);
    int d_out = OutDegree(node_i);
    vector<double> z_stars;
    z_stars.push_back(exp((double)-d_in) - exp(-((double) d_in+ (double) ind)));
    z_stars.push_back(exp((double)-d_out) - exp(-((double)d_out+(double)ind)));
    return z_stars;
}


double Two_Path_Change(std::string node_i, std::string node_j,int xij)
{
    double z_two_path = 0;
    double x_p;
    if (xij == 0) x_p = -1.0;
    else x_p = 1.0;

    std::string Mode_Out = "OUT";
    std::string Mode_In = "IN";

    std::set<std::string> out_neighbors_j = Neighbors(node_j,Mode_Out);
    std::set<std::string> in_neighbors_i = Neighbors(node_i,Mode_In);
    int path_i = in_neighbors_i.size();
    int path_j = out_neighbors_j.size();
    z_two_path += x_p*2.0*(1 - pow(0.5, (double)path_i)) + x_p*2.0*(1 - pow(0.5, (double)path_j));


    return z_two_path;
}

// ============================================== End of functions to calculate the change statistics ============================================================================

// =================================================================== Fixed density MCMC =========================================================================================

std::pair<VectorXd, std::unordered_map <std::string, struct Link> > IFD_Sampling(vector<double> current_configuration)
{
    int N_A, N_D,xij;
    bool test_A, test_D;
	int z_edge = 0;
	int z_reciprocity = 0;
	int z_sector = 0;
	double z_geographic = 0.0;
	int z_bank = 0;
	int z_employee = 0;

	vector<double> z_profit;
    z_profit.push_back(0.0);
    z_profit.push_back(0.0);
    z_profit.push_back(0.0);
	vector<double> z_triangle;
    z_triangle.push_back(0.0);
    z_triangle.push_back(0.0);
    z_triangle.push_back(0.0);
    z_triangle.push_back(0.0);
    z_triangle.push_back(0.0);
	vector<double> z_stars;
    z_stars.push_back(0.0);
    z_stars.push_back(0.0);
	double z_path = 0.0;

    std::string source_node, target_node;
    for (int i = 0; i < IFD_size; ++i)
    {
        test_A = false;
        while ((test_A == false))
	    {
            bool no_link = false;
            while (no_link == false)
            {
        		int random_source = rand() % network.size();
        		int random_target = rand() % network.size();
        		while (random_source == random_target) random_target = rand() % network.size();

        		std::unordered_map<std::string, struct Link>::iterator itr_source_node = network.begin();
        		std::unordered_map<std::string, struct Link>::iterator itr_target_node = network.begin();
                std::advance(itr_source_node, random_source);
        		std::advance(itr_target_node, random_target);

    			source_node = (*itr_source_node).first;
    			target_node = (*itr_target_node).first;

        		if ((Adjacency(source_node, target_node) == 0) && ((source_node != "") && (target_node != ""))) no_link = true;
            }

            xij = 1;
			network[source_node].out_neighbors.insert(target_node);
			network[target_node].in_neighbors.insert(source_node);

	        z_edge = Edge_Change(xij);
	        z_reciprocity = Reciprocity_Change(source_node, target_node,xij);
	        z_sector = Sector_Change(source_node, target_node,xij);
	        z_bank = Bank_Change( source_node, target_node,xij);
	        z_employee = Employee_Change(source_node, target_node,xij);
	        z_geographic = Geographic_Change(source_node, target_node,xij);
	        z_profit = Profit_Change(source_node, target_node,xij) ;
	        z_triangle = Triangles_Change(source_node, target_node,xij) ;
	        z_stars = Stars_Change(source_node, target_node,xij);
	        z_path = z_triangle[4];

	        double P_1 = exp(theta(0)*z_edge + z_reciprocity*theta(1)+ z_sector*theta(2)+ z_geographic*theta(3)+ z_profit[0]*theta(4)+ z_profit[1]*theta(5)
	       		+ z_profit[2]*theta(6)+z_stars[0]*theta(7) + z_stars[1]*theta(8) + z_path*theta(9)+ z_triangle[0]*theta(10)+ z_triangle[1]*theta(11)
	       		+ z_triangle[2]*theta(12)+ z_triangle[3]*theta(13) + z_bank*theta(14) + z_employee*theta(15));

	        double P = (double) rand() / (RAND_MAX);
	        double Hasting_ratio = min(1.0,P_1);

	        if (P<Hasting_ratio)
	        {
	        	vector<double> change_statistics;
	        	change_statistics.push_back(z_edge);
                change_statistics.push_back(z_reciprocity);
                change_statistics.push_back(z_sector);
                change_statistics.push_back(z_geographic);
                change_statistics.push_back(z_profit[0]);
                change_statistics.push_back(z_profit[1]);
                change_statistics.push_back(z_profit[2]);
                change_statistics.push_back(z_stars[0]);
                change_statistics.push_back(z_stars[1]);
                change_statistics.push_back(z_path);
                change_statistics.push_back(z_triangle[0]);
                change_statistics.push_back(z_triangle[1]);
                change_statistics.push_back(z_triangle[2]);
                change_statistics.push_back(z_triangle[3]);
                change_statistics.push_back(z_bank);
                change_statistics.push_back(z_employee);

                std::transform (current_configuration.begin(), current_configuration.end(), change_statistics.begin(), current_configuration.begin(), std::plus<double>());
	        	test_A = true;
	        }

	        else
	        {
				network[source_node].out_neighbors.erase(target_node);
				network[target_node].in_neighbors.erase(source_node);
	        }
	    }

        test_D = false;
        while ((test_D == false))
        {
        	bool out_neighbor_test = false;
        	bool no_link = false;
        	std::set<std::string> out_neighbor;
        	while (no_link == false)
        	{
            	while(out_neighbor_test == false)
            	{
                	int random_source = rand() % network.size();
                    std::unordered_map<std::string, struct Link>::iterator itr_source_node = network.begin();
                	std::advance(itr_source_node, random_source);
                	source_node = (*itr_source_node).first;
                	out_neighbor  = network[source_node].out_neighbors;
                	if (out_neighbor.size() != 0) out_neighbor_test = true;
            	}

            	int random_target = rand() % out_neighbor.size();
                std::set<std::string>::iterator itr_target_node = out_neighbor.begin();
            	std::advance(itr_target_node, random_target);
            	target_node = *itr_target_node;

        		if ((Adjacency(source_node, target_node) == 1) && ((source_node != "") && (target_node != ""))) no_link = true;
        	}


        	xij = 0;

			network[source_node].out_neighbors.erase(target_node);
			network[target_node].in_neighbors.erase(source_node);

	        z_edge = Edge_Change(xij);
	        z_reciprocity = Reciprocity_Change(source_node, target_node,xij);
	        z_sector = Sector_Change(source_node, target_node,xij);
	        z_bank = Bank_Change(source_node, target_node,xij);
	        z_employee = Employee_Change(source_node, target_node,xij);
	        z_geographic = Geographic_Change(source_node, target_node,xij);
	        z_profit = Profit_Change(source_node, target_node,xij) ;
	        z_triangle = Triangles_Change(source_node, target_node,xij) ;
	        z_stars = Stars_Change(source_node, target_node,xij);
	        z_path = z_triangle[4];

	        double P_1 = exp(theta(0)*z_edge + z_reciprocity*theta(1)+ z_sector*theta(2)+ z_geographic*theta(3)+ z_profit[0]*theta(4)+ z_profit[1]*theta(5)
	        		+ z_profit[2]*theta(6)+z_stars[0]*theta(7) + z_stars[1]*theta(8) + z_path*theta(9)+ z_triangle[0]*theta(10)+ z_triangle[1]*theta(11)
	        		+ z_triangle[2]*theta(12)+ z_triangle[3]*theta(13) + z_bank*theta(14) + z_employee*theta(15));
	        double P = (double) rand() / (RAND_MAX);
	        double Hasting_ratio = min(1.0,P_1);

	        if (P<Hasting_ratio)
	        {
                vector<double> change_statistics;
	        	change_statistics.push_back(z_edge);
	        	change_statistics.push_back(z_reciprocity);
                change_statistics.push_back(z_sector);
                change_statistics.push_back(z_geographic);
                change_statistics.push_back(z_profit[0]);
                change_statistics.push_back(z_profit[1]);
                change_statistics.push_back(z_profit[2]);
                change_statistics.push_back(z_stars[0]);
                change_statistics.push_back(z_stars[1]);
                change_statistics.push_back(z_path);
                change_statistics.push_back(z_triangle[0]);
                change_statistics.push_back(z_triangle[1]);
                change_statistics.push_back(z_triangle[2]);
                change_statistics.push_back(z_triangle[3]);
                change_statistics.push_back(z_bank);
                change_statistics.push_back(z_employee);

                std::transform (current_configuration.begin(), current_configuration.end(), change_statistics.begin(), current_configuration.begin(), std::plus<double>());
	        	test_D = true;
	        }

	        else
	        {
    			network[source_node].out_neighbors.insert(target_node);
    			network[target_node].in_neighbors.insert(source_node);
	        }
        }
    }
	VectorXd Z_s(number_variables);
	for (int i = 0 ; i < number_variables; ++i)
	Z_s(i) = current_configuration[i];
	return std::make_pair(Z_s, network);
}

// =================================================================== End of Fixed density MCMC =========================================================================================

// =================================================================== Arithmetics operations =========================================================================================

double sum(vector<double> a)
{
	double s = 0;
	for (int i = 0; i < a.size(); i++)
	{
		s += a[i];
	}
	return s;
}

double mean(vector<double> a)
{
	return sum(a) / a.size();
}

double sqsum(vector<double> a)
{
	double s = 0;
	for (int i = 0; i < a.size(); i++)
	{
		s += pow(a[i], 2);
	}
	return s;
}

double stdev(vector<double> nums)
{
	double N = nums.size();
	return pow(sqsum(nums) / N - pow(sum(nums) / N, 2), 0.5);
}

// =================================================================== End of Arithmetics operations =========================================================================================

// =================================================================== Robins Monro Stochastic Approximation =============================================================================


VectorXd MLE_Robins_Monro(vector<double>& configuration)
{
	VectorXd Z_obs(number_variables);
    std::unordered_map <std::string, struct Link> sampled_network;
    VectorXd Z_s;
    std::pair<VectorXd, std::unordered_map <std::string, struct Link> > Simulation;
	for (int i = 0 ; i < number_variables; ++i)
	Z_obs(i) = configuration[i];
	VectorXd E(number_variables);
	E << 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	MatrixXd M(initial_simulation, number_variables);
	MatrixXd D(number_variables,number_variables);
	for (int k = 0; k <initial_simulation; ++k)
	{
		Simulation = IFD_Sampling(configuration);
        Z_s = Simulation.first;
        sampled_network = Simulation.second;
		M.row(k) = Z_s;
		for (int i = 0; i<number_variables;++i) E[i] += Z_s[i];
	}
	double weight_parameters = 1.0/initial_simulation;
    E = weight_parameters*E;
    D = weight_parameters*M.transpose()*M - E*E.transpose();

    FullPivLU< MatrixXd > D_lu(D);
	MatrixXd inv_D = D_lu.inverse();
    MatrixXd inv_D0(number_variables,number_variables);
    for (int i = 0; i <number_variables; ++i)
    {
    	for (int j = 0; j <number_variables; ++j)
    	{
    		if (i==j)
    		{
    			if (D.row(i)[j] != 0) inv_D0.row(i)[j] = 1/D.row(i)[j];
    			else inv_D0.row(i)[j] = 0;
    		}
    		else inv_D0.row(i)[j] = 0;
    	}
    }

    theta -= a_r*inv_D0*(E - Z_obs);

    int nbr_sub_phases;

    typedef Matrix<double,16,1> vec_theta;
    std::vector<vec_theta> Theta_Collection;
    VectorXd theta_average(number_variables);

	for (int m = 0; m<burnin; ++m)
	{
		Simulation = IFD_Sampling(configuration);
        Z_s = Simulation.first;
        sampled_network = Simulation.second;
		theta -= a_r*inv_D0*(Z_s - Z_obs);
	}

    int r = 0;
    bool convergence_test = false;
    int countconv = 0;
    while (r<=number_phases && convergence_test == false)
    {
    	nbr_sub_phases = 10;
    	theta_average << 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
        vector<vector<double> > Statistics_Collection(number_variables, vector<double>(1));
    	for (int m = 0; m<nbr_sub_phases; ++m)
    	{
    		Simulation = IFD_Sampling(configuration);
            Z_s = Simulation.first;
            sampled_network = Simulation.second;
    		theta -= a_r*inv_D0*(Z_s - Z_obs);
    		theta_average+=theta;
    		Theta_Collection.push_back(theta);
    		for (int i = 0; i<number_variables; ++i) Statistics_Collection[i].push_back(Z_s[i]);
    	}

    	a_r = 0.5*a_r;
    	theta = theta_average/nbr_sub_phases;

		if (r > 3)
		{
			int n_conv = 0;
			for (int i= 0; i<number_variables; ++i)
			{
				if (abs((mean(Statistics_Collection[i]) - Z_obs[i])/stdev(Statistics_Collection[i])) <= 2)
					n_conv++;
			}
			if (n_conv == number_variables) convergence_test = true;
		}
    	++r;
    }

	ofstream theta_data;
    std::string path = "Theta_Collection.txt";
	theta_data.open(path);
	for (int i=0; i<Theta_Collection.size(); ++i)
	{
		for (int j = 0; j < Theta_Collection[i].size(); ++j)
		{
		    theta_data << Theta_Collection[i][j] << '\t' ;
		}
		theta_data << '\n' ;
	}
	theta_data.close();
	return theta;
}

// =================================================================== End of Robins Monro Stochastic Approximation ==========================================================================

// ==================================================== Simulate network based on the estimated Theta -- To use for GoF analysis =============================================================


void Save_Network(std::unordered_map<std::string, struct Link>& SimNetwork, int& i_path)
{
    ofstream network_data;
    std::string path1 = "network_sim";
    std::string path2 = ".txt";
    std::string path = path1 + std::to_string(i_path) + path2;
    network_data.open(path);

    for (std::unordered_map<std::string, struct Link>::iterator source_node=SimNetwork.begin(); source_node!=SimNetwork.end();++source_node)
    {
    	for (std::set<std::string>::iterator target_node=SimNetwork[(*source_node).first].out_neighbors.begin(); target_node!=SimNetwork[(*source_node).first].out_neighbors.end();++target_node)
    	{
    		network_data << (*source_node).first  << '\t' ;
    		network_data << *target_node  << '\t' ;
    		network_data << '\n' ;
    	}
    }

    network_data.close();
}

// ============================================================== Simulate network based on the estimated Theta -- End ======================================================================

int main() {
	theta << 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	Upload_Data();
	vector<double> current_configuration = Network_Configuration();
	VectorXd result = MLE_Robins_Monro(current_configuration);

    std::pair<VectorXd, std::unordered_map <std::string, struct Link> > Simulation;
    for (int i = 0; i < NumNetworks; ++i)
    {
        Simulation = IFD_Sampling(current_configuration);
        Save_Network(Simulation.second,i);
    }

	return 0;
}
