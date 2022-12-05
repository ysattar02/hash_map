#ifndef DNADB_H
#define DNADB_H
#include <iostream>
#include <string>
#include "math.h"
using namespace std;
class Grader;   // forward declaration, will be used for grdaing
class Tester;   // forward declaration, will be used for testing
class DNA;      // forward declaration
class DnaDb;    // forward declaration
const int MINLOCID = 1000;
const int MAXLOCID = 9999;
const int MINPRIME = 101;   // Min size for hash table
const int MAXPRIME = 99991; // Max size for hash table
#define EMPTY DNA("")
#define DELETED DNA("DELETED")
#define DELETEDKEY "DELETED"
typedef unsigned int (*hash_fn)(string); // declaration of hash function
const int MAX = 4;
const char ALPHA[MAX] = {'A', 'C', 'G', 'T'};

class DNA{
public:
    friend class Grader;
    friend class Tester;
    friend class DnaDb;
    DNA(string sequence="", int location=0); // Constructor
    string getSequence() const;              // Returns the key
    int getLocId() const;
    // Overloaded assignment operator
    const DNA& operator=(const DNA& rhs);
    // Overloaded insertion operator
    friend ostream& operator<<(ostream& sout, const DNA &dna );
    // Overloaded equality operator
    friend bool operator==(const DNA& lhs, const DNA& rhs);
private:
    string m_sequence;  // this is the object key
    int m_location;     // some info
};

class DnaDb{
public:
    friend class Grader;
    friend class Tester;
    DnaDb(int size, hash_fn hash);
    ~DnaDb();
    // Returns Load factor of the new table
    float lambda() const;
    // Returns the ratio of deleted slots in the new table
    float deletedRatio() const;
    // insert only happens in the new table
    bool insert(DNA dna);
    // remove can happen from either table
    bool remove(DNA dna);
    // find can happen in either table
    DNA getDNA(string sequence, int location);
    void dump() const;

private:
    hash_fn         m_hash;         // hash function

    DNA*            m_currentTable; // hash table
    unsigned int    m_currentCap;   // hash table size
    unsigned int    m_currentSize;  // current number of entries
    // m_currentSize includes deleted entries
    unsigned int    m_currNumDeleted;// number of deleted entries

    DNA*            m_oldTable;     // hash table
    unsigned int    m_oldCap;       // hash table size
    unsigned int    m_oldSize;      // current number of entries
    // m_oldSize includes deleted entries
    unsigned int    m_oldNumDeleted;// number of deleted entries

    //private helper functions
    bool isPrime(int number);
    int findNextPrime(int current);

    /******************************************
    * Private function declarations go here! *
    ******************************************/
    enum class REHASH_STATUS { NOT_REHASHING, QUARTER, HALF, THREE_QUARTER };
    REHASH_STATUS rehash_status;
    unsigned int get_index_cur(DNA dna, bool deleted_empty) const;
    unsigned int get_index_old(DNA dna, bool deleted_empty) const;
    void rehash();
    friend class Tester;
};
#endif