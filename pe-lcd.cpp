#include <iostream>
#include <ctime>
#include <sstream>
#include <cstdlib>
#include<cstdlib>
#include <fstream>
#include <chrono>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <iterator>
#include <algorithm>
#include <math.h>
#include <iomanip>
#include <ios>
#include <string>
#include <set>
#include <map>
#include <list>
#include <numeric>
using namespace std;
class Graph
    {
		long ecount = 0;
    public:
        int max_degree=0;
        set<int> V;//only vertices with positive degree are stored
        map<int, set<int> > neighbours;
        map<int,int>Centnode;
        map<int,int>Centnode1;
        void read_edgelist(string&);
        inline int vcount(){ return V.size(); }
        inline long get_ecount() { return ecount; }
        void print_graph();
    };

    void Graph::read_edgelist(string& edgefile)
 {
    ifstream fin;
    fin.open(edgefile);
    if(!fin.is_open())
    {
        cout<<"The file containing the edgelist could not be opened."<<endl;
        exit(1);
    }

    string line;
    while ( getline(fin, line ) )
    {
        istringstream is( line );
        int u, v;
        is>>u;
        is>>v;
        if(u != v and neighbours[u].find(v) == neighbours[u].end())
        {
            neighbours[u].insert(v);
            neighbours[v].insert(u);
            ecount++;
        }
        V.insert(u);
        V.insert(v);
    }

    set<int>::iterator si,sj,st;
    for(st=V.begin();st!=V.end();++st)
    {
        if(neighbours[*st].size()>=1){
             int i=*st;
             int node=0;
             for(si=neighbours[i].begin();si!=neighbours[i].end();++si){
                   for(sj=neighbours[i].begin();sj!=neighbours[i].end();++sj){
                   if(si!=sj&&neighbours[*si].find(*sj)!=neighbours[*si].end())
                node++;
               }
           }
       Centnode.insert({i,neighbours[i].size()+(node/2)});
      }
    }

    Centnode1=Centnode;
}

void Graph::print_graph()
{

    cout<<endl<<"--------"<<endl;
    cout<<"No. of vertices: "<<vcount()<<endl;
    cout<<"No. of edges: "<<ecount<<endl;
    cout<<"\n"<<"vertex"<<setw(8)<<"   neighbours"<<endl;
    map<int, set<int> >::iterator it;
  /* for(it=neighbours.begin(); it != neighbours.end(); ++it)
      {
        cout<<it->first<<setw(5)<<"--> ";
        set<int>::iterator si;
        for(si = it->second.begin(); si != it->second.end(); ++si)
            cout<<*si<<" ";
        cout<<endl;
    }*/




}
class community
{
	int interior_ecount = 0;
	public :
	set<int> interior;
//public:
    int unstability ;
    int extdeg = 0;
	set<int> boundary;
	set<int>neighbors;
	bool stable = false;
	set<int> nearest_coms;
	inline int int_vcount(){ return (interior.size()); }
	inline int vcount() { return (interior.size()+boundary.size()); }
	inline int int_ecount() { return interior_ecount; }
	void get_new_nodes(Graph&, const set<int>&, set<int>&, float);
	void update_interior(Graph&, set<int>& );
	void update(Graph&, community&, int );
	void set_stability(Graph&, set<int>& );
	float overlap(Graph&, community& );
	void print();
	void write(ostream& );	//can be used for printing to the terminal and to file
};
//other function declarations
void usage(void);
float interaction_coefficient(Graph&, set<int>&, int);
float interaction_coefficient(Graph&, set<int>&, int, int);
float overlap(Graph&, community&, community&);
void set_nearest_coms(Graph&, map<int, community>&);
void merge_communities(Graph&, map<int, community>& );
void write_communities(string&, map<int, community>&);
void print_communities(map<int, community>&, set<int>&);
int rand(int, int);
void print_set(const set<int>&);
void print_map(map<int, int>&);

//community methods definitions
void community::update_interior(Graph& g, set<int>& nbrs_c)
{
    set<int>::iterator si, sj;
    //cout<<"all boundary points: "; print_set(boundary);
    si = boundary.begin();
    //cout<<"curr boundary point = "<<*si<<endl;
    do
    {
        bool all_nbrs_in_c = true;
        for(sj = g.neighbours[*si].begin(); sj != g.neighbours[*si].end(); ++sj)
            if(nbrs_c.find(*sj) != nbrs_c.end())
            {
                all_nbrs_in_c = false;
                break;
            }
        //cout<<"all nbrs in c = "<<all_nbrs_in_c<<endl;
        if(all_nbrs_in_c)
        {
            for(sj = g.neighbours[*si].begin(); sj != g.neighbours[*si].end(); ++sj)
                if(interior.find(*sj) != interior.end())
                    interior_ecount++;
            //cout<<"interior_ecount = "<<interior_ecount<<endl;
            interior.insert(*si);
            si = boundary.erase(si);   //remove current and point to next
        }
        else
            ++si;
    }while(si != boundary.end());
}



void community::write(ostream& location)
{
    interior.insert(boundary.begin(), boundary.end());
    copy(interior.begin(), interior.end(), ostream_iterator<int>(location, " "));
}
float interaction_coefficient(Graph& g, set<int>& s, int v)
{
   float icoeff = 0;
   set<int>::iterator si;
    for(si = g.neighbours[v].begin(); si != g.neighbours[v].end(); ++si)
        if(s.find(*si) != s.end())
            icoeff += 1;
    icoeff /= g.neighbours[v].size();
    return icoeff;
}
float interaction_coefficient(Graph& g, set<int>& s, int u, int v)
{
    float icoeff = 0, indeg_u = 0, indeg_v = 0;
    set<int>::iterator si;
    for(si = g.neighbours[u].begin(); si != g.neighbours[u].end(); ++si)
        if(s.find(*si) != s.end())
            indeg_u += 1;
    for(si = g.neighbours[v].begin(); si != g.neighbours[v].end(); ++si)
        if(s.find(*si) != s.end())
            indeg_v += 1;
    icoeff = 1 + min(indeg_u, indeg_v);
    icoeff /= max(g.neighbours[u].size(), g.neighbours[v].size()); //counting the edge uv only once
    return icoeff;
}

void community::get_new_nodes(Graph& g, const set<int>& nbrs_c, set<int>& new_nodes, float rho_0)
{
    set<int> curr_nbrs = nbrs_c;
    set<int>::iterator si, sj;
    float icoeff;
    while(!curr_nbrs.empty())
    {
        si = curr_nbrs.begin();
        //cout<<"curr nbr = "<<*si<< "'s icoeff = ";
        icoeff = interaction_coefficient(g, boundary, *si);
        //cout<<icoeff<<endl;
        if( icoeff >= rho_0)
            new_nodes.insert(*si);
        else
        {
            for(sj = g.neighbours[*si].begin(); sj != g.neighbours[*si].end(); ++sj)
                if(nbrs_c.find(*sj) != nbrs_c.end() && *sj > *si)
                {
                    //cout<<"edge "<<*si<<"-"<<*sj<<"'s icoeff = ";
                    icoeff = interaction_coefficient(g, boundary, *si, *sj);
                    //cout<<icoeff<<endl;
                    if(icoeff >= rho_0)
                    {
                        new_nodes.insert(*si);
                        new_nodes.insert(*sj);
                        curr_nbrs.erase(*sj);
                    }
                }
        }
        curr_nbrs.erase(*si);
    }
    if(new_nodes.empty() and interior.empty())
    {
        set<int> c_exterior = boundary;
        c_exterior.insert(nbrs_c.begin(), nbrs_c.end());
        //cout<<"c's exterior "<<endl; print_set(c_exterior);
        for(si = nbrs_c.begin(); si != nbrs_c.end(); ++si)
            if(interaction_coefficient(g, c_exterior, *si) >= 0.3)
                new_nodes.insert(*si);
    }
}




void write_communities(string& comfile, map<int, community>& Cover)
{
    map<int, community>::iterator msi;
    set<int>::iterator si;
    ofstream fout(comfile);
    if(!fout.is_open())
    {
        cout<<"Destination file for communities could not be opened."<<endl;
        exit(1);
    }
    for(msi = Cover.begin(); msi != Cover.end(); ++msi)
    {
        msi->second.write(fout);
        fout<<endl;
    }
}

void print_set(const set<int>& s)
{
    set<int>::iterator sitr1;
    for(sitr1 = s.begin(); sitr1 != s.end(); ++sitr1)
            cout<<*sitr1<<" ";
    cout<<endl;
}


pair<int, int>LargestValue(map<int, int> sampleMap)
{
	pair<int, int> MaxValue = *sampleMap.begin();
	map<int, int>::iterator currentEntry;
	for (currentEntry = sampleMap.begin();currentEntry != sampleMap.end();++currentEntry) {
		if (currentEntry->second> MaxValue.second) {
             MaxValue= make_pair(currentEntry->first,currentEntry->second);
		   }
	    }


    return MaxValue;
}




int main(int argc, char* argv[])
{
    std::chrono::time_point<std::chrono::system_clock> start_time, end_time;
    start_time = std::chrono::system_clock::now();  //time starts
    float  rho_0 = 0.30;
    Graph g;
    string network;

     if(argc < 2 || argc > 4)
     {
	cout<<"Please see README file."<<endl;
	exit(1);
     }
    if(argv[1][0] == '-')
    {
	cout<<"Please see README file."<<endl;
	exit(1);
     }
    else
    {
        network = string(argv[1]);
    }
     int k=2;
    while(k <= argc-1)
    {
        string arg = string(argv[k]);
        if(arg == "-rh")
              {
               istringstream is(argv[k+1]);
               is>>rho_0;
               if( rho_0 < 0 || rho_0 >= 1)
		  {
		    cout<<"Please see README file."<<endl;
		    exit(1);
		   }
		  k+=2;
	       }
	    else
              {
		    cout<<"PLEASE see README file."<<endl;
		    exit(1);
		   }
      }

    g.read_edgelist(network);
    cout<<"-------------------------------------------------------------------"<<endl;
    cout<<"Pure Expansion based Local Community Detection (PE-LCD) "<<endl;
    cout<<"Authors: Abhinav Kumar, Pawan Kumar and Ravins Dohare"<<endl;
    cout<<"Email: abhinavkumar080395@gmail.com, pkumariitd@gmail.com, ravinsdohare@gmail.com"<<endl;
    cout<<"-------------------------------------------------------------------"<<endl;
    map<int, community> Cover;
    map<int, set<int> > members;
    set<int> seeds;
    set<int>::iterator si, sj;
    map<int, community>::iterator mis;
    map<int,set<int>>::iterator miis;
    //cout<<"Initial network size  ="<<g.Centnode.size()<<endl;
   int i = 1;
    while(!g.Centnode.empty())
    {
        pair<int,int>yup=LargestValue(g.Centnode);
       int r=yup.first;
      // cout<<endl<<"seed node ="<<r<<endl;
        community c;
        c.boundary.insert(r);
        c.extdeg = g.neighbours[r].size();
        g.Centnode.erase(r);
        set<int> nbrs_c = g.neighbours[r];
        set<int> new_nodes;
        do
        {
            new_nodes.clear();
            c.get_new_nodes(g, nbrs_c, new_nodes, rho_0);
            //print_set(new_nodes);

            for(si = new_nodes.begin(); si != new_nodes.end(); ++si)
            {
                nbrs_c.erase(*si);
                for(sj = g.neighbours[*si].begin(); sj != g.neighbours[*si].end(); ++sj)
                    if(c.boundary.find(*sj) != c.boundary.end())
                        c.extdeg--;
                    else if(new_nodes.find(*sj) == new_nodes.end())
                    {
                        c.extdeg++;
                        nbrs_c.insert(*sj);
                    }
            }
            c.boundary.insert(new_nodes.begin(), new_nodes.end());

            if(!new_nodes.empty())
                c.update_interior(g, nbrs_c);
            for(si = new_nodes.begin(); si != new_nodes.end(); ++si)
                g.Centnode.erase(*si);
        }while(!new_nodes.empty());
       for (si=nbrs_c.begin();si!=nbrs_c.end();++si)
            c.neighbors.insert(*si);
            Cover[i] = c;
              ++i;
        //cout<<"Remaining network size = "<<g.Centnode.size()<<endl;
    }
//cout<<"Writing final communities "<<endl;

ostringstream comfile;
    comfile<<"./pelcd-coms.txt";
    string c_file = comfile.str();
    	write_communities(c_file, Cover);
 end_time = std::chrono::system_clock::now();  //time ends
    std::chrono::duration<double> elapsed_seconds = end_time-start_time;
    cout.setf(ios::left, ios::adjustfield);
    cout<<"-------------------------------------------------------------------"<<endl;
    cout<<setw(33)<<"Network file"<<"= "<<setw(20)<<network<<endl;
        cout<<setw(33)<<"No. of vertices"<<"= "<<setw(20)<<g.vcount()<<endl;
        cout<<setw(33)<<"No. of edges"<<"= "<<setw(20)<<g.get_ecount()<<endl;
        cout<<setw(33)<<"Total communities"<<"= "<<Cover.size()<<endl;
        cout<<setw(33)<<"Community file"<<"= "<<setw(20)<<comfile.str()<<endl;
        cout<<setw(33)<<"Time elapsed"<<"= "<<elapsed_seconds.count()<<"s\n";
        cout<<"-------------------------------------------------------------------"<<endl;
    return 0;
}
