#include  <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <list>
#include <ctime>
#include <algorithm>
#include <math.h>
#include <stdint.h>

using namespace std;

#define	       PH_SIZE 	(1<<22)
#define     WH_SIZE         10000000
int kSize=11;
int winSize=2048;
int l_ref;
string s_ref;
pair <int,int> ph[PH_SIZE];
pair<uint32_t, int> wh[WH_SIZE];
int windows[10000];


/*to read reference file and
	 build the hsh map*/
//convert kmer to bit
/*uint32_t toBit(string s)
{
    uint32_t bitnum = 0, shift;
    for (size_t i = 0; i < kSize; i++)
    {
        shift = 2* kSize-2-2*(i%( kSize));
        switch (s[i])
        {
        case 'A':
            bitnum |= 0x0 << shift;
            break;
        case 'C':
            bitnum |= 0x1 << shift;
            break;
        case 'G':
            bitnum |= 0x2 << shift;
            break;
        case 'T':
            bitnum |= 0x3 << shift;
            break;
        default:
            cout << " invalid base\n";
            return 0;
        }
    }
    return bitnum;
}*/

uint32_t toBit(string s)
{
    uint32_t num = 0;
    for (int i = 0; i < s.size(); i++)
    {
        shift = 2* kSize-2-2*(i%( kSize));
        switch (s[i])
        {
        case 'A':
            num = (num<<2)+0;
            break;
        case 'C':
            num = (num<<2)+1;
            break;
        case 'G':
            num = (num<<2)+2;
            break;
        case 'T':
            num = (num<<2)+3;
            break;
        default:
            cout << " invalide char\n";
            return 0;
        }
    }
    return num;
}

bool cmp(pair<uint32_t, int> a,pair<uint32_t, int> b)
{
    return a.first<b.first;
}
//build hash map
void createHash()
{

    //read reference file
    



        //define a  hash table
        int i,j,len=0;
        for(i=0; i<l_ref-kSize; i++)
        {
            if(i<1024)
            {
                wh[len++]=make_pair(toBit(s_ref.substr(i,kSize)),1);
            }
            else
            {
                j=i/1024;
                wh[len++]=make_pair(toBit(s_ref.substr(i,kSize)),j);
                wh[len++]=make_pair(toBit(s_ref.substr(i,kSize)),j+1);
            }
        }

        sort(wh,wh+len,cmp);


       /* int i,j,k,n=l_ref/1024,len=0;
        for(i=0; i<1024; i++)
        {
            wh[len++]=make_pair(toBit(s_ref.substr(i,kSize)),1);
        }
        for(j=1; j<n; j++)
        {
            k=i;
            for(i=k; i<k+1024; i++)
            {
                wh[len++]=make_pair(toBit(s_ref.substr(i,kSize)),j);
                wh[len++]=make_pair(toBit(s_ref.substr(i,kSize)),j+1);
            }
        }
        for(; i<l_ref-kSize; i++)
        {
            wh[len++]=make_pair(toBit(s_ref.substr(i,kSize)),j);
        }*/

        //sort(wh,wh+len,cmp);


          /*   time_t start,stop;
        start = time(NULL);
         stop = time(NULL);
        printf("Use Time:%ld\n",(stop-start));*/
        i=0;
        while(i<len)
        {
            j=wh[i].first;
            ph[j].first=i;
            while(wh[i].first==j && i<len)
            {
                i++;
            }
            ph[j].second=i;
        }
}


//find the window which has the most kmers
int findWindow(string s)
{
    int winNum,times,bitNum=0;
    int i,j,slen=s.length();
    for(i=0; i<slen-kSize; i++)
    {
        bitNum=toBit(s.substr(i,kSize));
        for( j=ph[bitNum].first; j<ph[bitNum].second; j++)
        {
            windows[wh[j].second]+=1;
        }
    }
    times=windows[0];
    winNum=0;
    for(i=0; i<l_ref/1024 ; i++)
    {
        if(times<windows[i])
        {
            times=windows[i];
            winNum=i;
        }
    }
    cout<<winNum<<endl;
    return winNum;
}

//alignment with read and window which is finded


int main()
{
    time_t start,stop;
    start = time(NULL);

    string line;
    ifstream in_ref("E.coli.fa");
    if(!in_ref)
    {
        cout<<"fail to open the file!"<<endl;
    }
    else
    {
        getline(in_ref,line);
        while(getline(in_ref,line))
        {
            s_ref+=line;
        }
        l_ref=s_ref.length();
        createHash();
        in_ref.close();
    }

    //creaste hash map
    

    //read reads file
    string s_read,line,s_read_mid;
    int l_read,readNum=0;
    ifstream in_read("E.coli-sim.fastq.coli-sim");
    if(!in_read)
    {
        cout<<"fail to open the file!"<<endl;
        return -1;
    }
    else
    {
        while(!in_read.eof())
        {
            readNum++;
            getline(in_read,line);
            getline(in_read,s_read);
            l_read=s_read.length();
            s_read_mid=s_read.substr(l_read/2,1000);
            getline(in_read,line);
            getline(in_read,line);
            cout<<"read"<<readNum<<":";
            findWindow(s_read_mid);
        }
        in_read.close();
    }

    stop = time(NULL);
    printf("Use Time:%ld\n",(stop-start));
    return 0;
}



