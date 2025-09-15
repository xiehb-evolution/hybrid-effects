#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <sstream>

using namespace std;

using std::ios;
using std::cout;
using std::endl;

vector<string> split(const string& s, const string& delim, const bool keep_empty = true)
{
    vector<string> result;
    if (delim.empty()) {
        result.push_back(s);
        return result;
    }
    string::const_iterator substart = s.begin();
    string::const_iterator subend;
    while (true)
    {
        subend = search(substart, s.end(), delim.begin(), delim.end());
        string temp(substart, subend);
        if(keep_empty || !temp.empty())
        {
            result.push_back(temp);
        }
        if (subend == s.end())
        {
            break;
        }
        substart = subend + delim.size();
    }
    return result;
}

void SplitString(const string& str, const string& delimiters, vector<string> &elems, bool skip_empty=false)
{
    string::size_type pos, prev = 0;
    while((pos=str.find_first_of(delimiters,prev))!=string::npos)
    {
        if(pos>prev)
        {
            if(skip_empty && 1== pos - prev)
                break;
            elems.emplace_back(str,prev,pos-prev);
        }
        prev = pos + 1;
    }
    if(prev < str.size())
        elems.emplace_back(str, prev, str.size() - prev);
}

vector<string> LoadFileLinesIntoVector(string filename)
{
    ifstream infile;
    string data;
    vector<string> retdata;

    infile.open(filename.c_str());
    while (getline(infile,data))
    {
        if(data.length()==0)
        {
            continue;
        }
        retdata.push_back(data);
    }
    infile.close();

    return retdata;
}

vector<int> FindSamplesInVCFHeader(string header, vector<string> sampleid)
{
    int i,j;
    vector<int> retpos;
    vector<string> cols;

    SplitString(header,"\t",cols,false);
    //cols = split(header,"\t",true);
    for(i=0;i<sampleid.size();i++)
    {
        for(j=0;j<cols.size();j++)
        {
            if(cols[j] == sampleid[i])
            {
                retpos.push_back(j);
                break;
            }
        }
        if(j==cols.size())
            retpos.push_back(-1);
    }
    return(retpos);
}

vector<int> GetSampleGenotypeCodeWithDepth(vector<string> genotype,vector<int> samplepos,int mindepth,int maxdepth)
{
    int i;
    vector<int> retval;
    vector<string> fieldinfo;
    char base1,base2;
    int depth;
    int samplesize;
    samplesize = samplepos.size();

    for(i=0;i<samplesize;i++)
    {
        fieldinfo.clear();
        SplitString(genotype[samplepos[i]],":",fieldinfo,false);
        depth = atoi(fieldinfo[2].c_str());
        if( depth<mindepth || depth>maxdepth )
            continue;
        base1 = genotype[samplepos[i]].c_str()[0];
        base2 = genotype[samplepos[i]].c_str()[2];
        if( base1 =='0'&& base2 =='0')
        {
          retval.push_back(0);
          continue;
        }
        if( base1=='1' && base2=='1')
        {
          retval.push_back(2);
          continue;
        }
        if((base1=='0'&& base2=='1')||(base1=='1'&& base2=='0'))
        {
          retval.push_back(1);
          continue;
        }
        retval.push_back(-1);
    }
    return(retval);
}


vector<int> CountSampleGenotype(vector<int> genotype)
{
    int i;
    vector<int> retval(3);
    int samplesize;
    samplesize = genotype.size();
    
    for(i=0;i<samplesize;i++)
    {
      if(genotype[i]>=0)
        retval[genotype[i]]++;
    }
    return(retval);
}

float CountMiorAlleleFreq(vector<int> genotype1, vector<int> genotype2)
{
    float allsize;
    float allelecount;
    float allelefreq;
    allsize = 0.0;
    allelecount = 0.0;
    int i;
    int samplesize;
    samplesize = genotype1.size();
    
    for(i=0;i<samplesize;i++)
    {
      if(genotype1[i]>=0)
      {
        allsize = allsize+2;
        allelecount += genotype1[i];
      }
    }
    samplesize = genotype2.size();
    for(i=0;i<genotype2.size();i++)
    {
      if(genotype2[i]>=0)
      {
        allsize = allsize+2;
        allelecount += genotype2[i];
      }
    }
    allelefreq = allelecount/allsize;
    if(allelefreq>0.5)
        allelefreq = 1.0 - allelefreq;
    return(allelefreq);
}

int main(int argc,char *argv[])
{
  if(argc!=6)
  {
    cout << "Usage: "<<argv[0]<<" windowsize mindepth maxdepth popfile1 popfile2\n";
    return 0;
  }
  string popfile1, popfile2;
  vector<string> pop1_individuals,pop2_individuals;
  vector<int> pop1_individuals_columns,pop2_individuals_columns;
  int pop1size, pop2size;
  int mindepth,maxdepth;
  vector<int> pop1_genotype,pop2_genotype;

  string lastchr, thischr;
  long windowsize;
  long lastwindow, thiswindow;

  int snpflag;
  int i, j;

  string linedata;
  vector<string> data_columns;
  vector<int> genotype1,genotype2;
  vector<int> genotypecount1,genotypecount2;

  int DAFid;
  float allelefreq;

  int snps[10];
  int homogenotypecount1[10],hetgenotypecount1[10];
  int homogenotypecount2[10],hetgenotypecount2[10];
  
  windowsize = atol(argv[1]);
  mindepth = atoi(argv[2]);
  maxdepth = atoi(argv[3]);
  popfile1 = argv[4];
  popfile2 = argv[5];

  pop1_individuals = LoadFileLinesIntoVector(popfile1);
  pop2_individuals = LoadFileLinesIntoVector(popfile2);

  while( getline(cin,linedata) )
  {
      if(linedata.find("#CHROM",0)==0)
          break;
  }
  pop1_individuals_columns = FindSamplesInVCFHeader(linedata, pop1_individuals);
  pop2_individuals_columns = FindSamplesInVCFHeader(linedata, pop2_individuals);

  pop1size = pop1_individuals_columns.size();
  pop2size = pop2_individuals_columns.size();
  
  lastwindow=0;
  thiswindow=0;
  lastchr="";

  cout << "Chr" << "\t" << "Position" << "\tPop1size\tPop2size\t" << "SNPs" << "\t" << "HomoSites1" << "\t" << "HetSites1" << "\t" << "HomoSites2" << "\t" << "HetSites2" << "\tHetRatio1" << "\tHetRatio2" << endl;
  
  while( getline(cin,linedata) )
  {
    //snpflag = 0;
    //data_columns = split(linedata,"\t",true);
    data_columns.clear();
    SplitString(linedata,"\t",data_columns,false);
    thischr = data_columns[0];
    thiswindow = atol(data_columns[1].c_str())/windowsize;
    if(lastchr=="")
    {
        lastchr = thischr;
        lastwindow = thiswindow;
        for(i=0; i<10; i++)
        {
            snps[i] = 0;
            homogenotypecount1[i] = 0;
            hetgenotypecount1[i] = 0;
            homogenotypecount2[i] = 0;
            hetgenotypecount2[i] = 0;
        }
     }

     if( (thischr!=lastchr) || (lastwindow!=thiswindow) )
     {
         if(snps[0]>0)
         {
            cout << lastchr << "\t" << lastwindow*windowsize << "\t" << pop1size << "\t" << pop2size;
            for(i=0; i<10; i++)
                cout << "\t" << snps[i] << "\t" << homogenotypecount1[i] << "\t" << hetgenotypecount1[i] << "\t" << homogenotypecount2[i] << "\t" << hetgenotypecount2[i] << "\t" << ((float)hetgenotypecount1[i])/(homogenotypecount1[i]+hetgenotypecount1[i]) << "\t" << ((float)hetgenotypecount2[i])/(homogenotypecount2[i]+hetgenotypecount2[i]);
            cout << endl;
         }

         lastchr = thischr;
         lastwindow = thiswindow;
         for(i=0; i<10; i++)
         {
            snps[i] = 0;
            homogenotypecount1[i] = 0;
            hetgenotypecount1[i] = 0;
            homogenotypecount2[i] = 0;
            hetgenotypecount2[i] = 0;
         }
     }
      
     genotype1 = GetSampleGenotypeCodeWithDepth(data_columns,pop1_individuals_columns,mindepth,maxdepth);
     genotype2 = GetSampleGenotypeCodeWithDepth(data_columns,pop2_individuals_columns,mindepth,maxdepth);

     genotypecount1 = CountSampleGenotype(genotype1);
     genotypecount2 = CountSampleGenotype(genotype2);
     if ( genotypecount1[0]==0 && genotypecount1[1]==0 && genotypecount2[0]==0 && genotypecount2[1]==0 )
        continue;
     if ( genotypecount1[2]==0 && genotypecount1[1]==0 && genotypecount2[2]==0 && genotypecount2[1]==0 )
        continue;
     if( (genotypecount1[0]==0 && genotypecount1[1]==0 && genotypecount1[2]==0 ) || (genotypecount2[0]==0 && genotypecount2[1]==0 && genotypecount2[2]==0 ) )
        continue;

     allelefreq = CountMiorAlleleFreq(genotype1, genotype2);
     DAFid = static_cast<int>(allelefreq/0.05);
     if(DAFid==10) DAFid = 9;

     homogenotypecount1[DAFid] += genotypecount1[0] + genotypecount1[2];
     hetgenotypecount1[DAFid] += genotypecount1[1];
     homogenotypecount2[DAFid] += genotypecount2[0] + genotypecount2[2];
     hetgenotypecount2[DAFid] += genotypecount2[1];
      
     snps[DAFid] ++;

  }
  if(snps[0]>0)
  {
       cout << lastchr << "\t" << lastwindow*windowsize << "\t" << pop1size << "\t" << pop2size;
       for(i=0; i<10; i++)
         cout << "\t" << snps[i] << "\t" << homogenotypecount1[i] << "\t" << hetgenotypecount1[i] << "\t" << homogenotypecount2[i] << "\t" << hetgenotypecount2[i] << "\t" << ((float)hetgenotypecount1[i])/(homogenotypecount1[i]+hetgenotypecount1[i]) << "\t" << ((float)hetgenotypecount2[i])/(homogenotypecount2[i]+hetgenotypecount2[i]);
       cout << endl;
  }
  
  return 1;
}

