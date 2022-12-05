#include "dnadb.h"
DnaDb::DnaDb(int size, hash_fn hash)
        :m_hash(hash), m_currentTable(nullptr), m_currentCap(0), m_currentSize(0), m_currNumDeleted(0),
         m_oldTable(nullptr), m_oldCap(0), m_oldSize(0), m_oldNumDeleted(0),
         rehash_status(REHASH_STATUS::NOT_REHASHING)
{
    //done
    if (size < MINPRIME) {
        size = MINPRIME;
    }
    else if (size > MAXPRIME) {
        size = MAXPRIME;
    }
    else if (!isPrime(size)) { //if not a prime number
        size = findNextPrime(size);
    }
    m_currentTable = new DNA[size];
    m_currentCap = size;
}

DnaDb::~DnaDb() {
    //done
    if (m_currentTable != nullptr) {
        delete[] m_currentTable; //cuz an array
        m_currentTable = nullptr;
        m_currentCap = 0;
        m_currentSize = 0;
        m_currNumDeleted = 0;
    }
    if (m_oldTable != nullptr) {
        delete[] m_oldTable; //cuz an array
        m_oldTable = nullptr;
        m_oldCap = 0;
        m_oldSize = 0;
        m_oldNumDeleted = 0;
    }
    m_hash = nullptr;
}

bool DnaDb::insert(DNA dna) {
    //done
    if (dna.m_location < MINLOCID || dna.m_location > MAXLOCID) {
        //bad location, reject insert operation
        return false;
    }
    unsigned int index = get_index_cur(dna, true);
    if (m_currentTable[index] == dna) {
        //already have this DNA
        return false;
    }
    //else not duplicate
    m_currentTable[index] = dna;
    m_currentSize++;
    if (rehash_status == REHASH_STATUS::NOT_REHASHING) {
        if (lambda() > .5f) { //floating type of .5
            rehash();
        }
    }
    else {
        rehash();
    }
    return true;
}

bool DnaDb::remove(DNA dna) {
    //done
    if (dna.m_location < MINLOCID || dna.m_location > MAXLOCID) { //if bad location id
        //bad location, reject insert operation
        return false;
    }
    unsigned int index = get_index_cur(dna, false);
    if (m_currentTable[index] == dna) {
        m_currentTable[index] = DELETED;    // DNA is in current table
        m_currNumDeleted++;
    }
    else if (m_oldTable == nullptr) {
        return false;   // DNA not in any table
    }
    else {
        index = get_index_old(dna, false);
        if (m_oldTable[index] == dna) {
            m_oldTable[index] = DELETED;    //DNA is in old table
            m_oldNumDeleted++;
        }
        else {
            return false;   // DNA not in any table
        }
    }
    if (rehash_status == REHASH_STATUS::NOT_REHASHING) {
        if (deletedRatio() > .8f) { //floating type of .8
            rehash();
        }
    }
    else {
        rehash();
    }
    return true;
}

DNA DnaDb::getDNA(string sequence, int location) {
    //done
    if (location < MINLOCID || location > MAXLOCID) { //if bad location id
        //bad location, reject insert operation
        return EMPTY;
    }
    DNA target = DNA(sequence, location);
    unsigned int index = get_index_cur(target, false);
    if (m_currentTable[index] == target) {
        return m_currentTable[index];
    }
    if (m_oldTable == nullptr) {
        return EMPTY;   // DNA not in any table
    }
    index = get_index_old(target, false);
    if (m_oldTable[index] == target) {
        return m_oldTable[index];
    }
    return EMPTY;
}

float DnaDb::lambda() const {
    //done
    //returns load factor
    if (m_currentCap == 0) {
        return 1;
    }
    //casting to float
    return float(m_currentSize) / float(m_currentCap);
}

float DnaDb::deletedRatio() const {
    //done
    if (m_currentSize == 0) {
        //safety check
        return 1;
    }
    return float(m_currNumDeleted) / m_currentSize;
}

void DnaDb::dump() const {
    //done
    cout << "Dump for current table: " << endl;
    if (m_currentTable != nullptr)
        for (int i = 0; i < m_currentCap; i++) {
            cout << "[" << i << "] : " << m_currentTable[i] << endl;
        }
    cout << "Dump for old table: " << endl;
    if (m_oldTable != nullptr)
        for (int i = 0; i < m_oldCap; i++) {
            cout << "[" << i << "] : " << m_oldTable[i] << endl;
        }
}

bool DnaDb::isPrime(int number) {
    //done
    bool result = true;
    for (int i = 2; i <= number / 2; ++i) {
        if (number % i == 0) {
            result = false;
            break;
        }
    }
    return result;
}

int DnaDb::findNextPrime(int current) {
    //done
    //we always stay within the range [MINPRIME-MAXPRIME]
    //the smallest prime starts at MINPRIME
    if (current < MINPRIME) current = MINPRIME - 1;
    for (int i = current; i < MAXPRIME; i++) {
        for (int j = 2; j * j <= i; j++) {
            if (i % j == 0)
                break;
            else if (j + 1 > sqrt(i) && i != current) {
                return i;
            }
        }
    }
    //if a user tries to go over MAXPRIME
    return MAXPRIME;
}

DNA::DNA(string sequence, int location) {
    //done
    if ((location >= MINLOCID && location <= MAXLOCID) ||
        (location == 0 && sequence == "DELETED")) {
        // this is a normal or a DELETED object
        m_sequence = sequence;
        m_location = location;
    }
    else {
        // this is the empty object
        m_sequence = "";
        m_location = 0;
    }
}

string DNA::getSequence() const {
    //done
    return m_sequence;
}

int DNA::getLocId() const {
    //done
    return m_location;
}

// Overloaded assignment operator
const DNA& DNA::operator=(const DNA& rhs) {
    if (this != &rhs) {
        m_sequence = rhs.m_sequence;
        m_location = rhs.m_location;
    }
    return *this;
}

// Overloaded insertion operator.  Prints DNA's sequence (key),
// and the location ID. This is a friend function in DNA class.
ostream& operator<<(ostream& sout, const DNA& dna) {
    //done
    if (!dna.m_sequence.empty())
        sout << dna.m_sequence << " (Location ID " << dna.m_location << ")";
    else
        sout << "";
    return sout;
}

// Overloaded equality operator. This is a friend function in DNA class.
// To test inequality we may negate the results of this operator.
bool operator==(const DNA& lhs, const DNA& rhs) {
    //done
    return ((lhs.m_sequence == rhs.m_sequence) && (lhs.m_location == rhs.m_location));
}

unsigned int DnaDb::get_index_cur(DNA dna, bool deleted_empty) const {
    //done
    unsigned int index = m_hash(dna.getSequence()) % m_currentCap;
    unsigned int temp = 1;
    while (!(m_currentTable[index] == dna)) {   // as long as DNA object in table is not
        // what were searching for
        if (m_currentTable[index] == EMPTY) {
            break;
        }
        if (deleted_empty && m_currentTable[index] == DELETED) {
            break;
        }
        //Quadratic Probing
        index += temp;
        index %= m_currentCap;
        temp += 2;
        temp %= m_currentCap;
    }
    return index;
}

unsigned int DnaDb::get_index_old(DNA dna, bool deleted_empty) const {
    //done
    unsigned int index = m_hash(dna.getSequence()) % m_oldCap;
    unsigned int temp = 1;
    while (!(m_oldTable[index] == dna)) {   // as long as DNA object in table is not
        // what were searching for
        if (m_oldTable[index] == EMPTY) {
            break;
        }
        if (deleted_empty && m_oldTable[index] == DELETED) {
            break;
        }
        //Quadratic Probing
        index += temp;
        index %= m_oldCap;
        temp += 2;
        temp %= m_oldCap;
    }
    return index;
}

void DnaDb::rehash() {
    //done
    int datapoints;
    if (rehash_status == REHASH_STATUS::NOT_REHASHING) {
        m_oldTable = m_currentTable;
        m_oldCap = m_currentCap;
        m_oldSize = m_currentSize;
        m_currentSize = 0;
        m_oldNumDeleted = m_currNumDeleted;
        m_currNumDeleted = 0;
        m_currentCap = findNextPrime(4 * (m_oldSize - m_oldNumDeleted));
        m_currentTable = new DNA[m_currentCap];
        datapoints = (m_oldSize - m_oldNumDeleted + 3) / 4;
        rehash_status = REHASH_STATUS::QUARTER; //updating rehash status
    }
    else if (rehash_status == REHASH_STATUS::QUARTER) {
        datapoints = (m_oldSize - m_oldNumDeleted + 2) / 3;
        rehash_status = REHASH_STATUS::HALF; //updating rehash status
    }
    else if (rehash_status == REHASH_STATUS::HALF) {
        datapoints = (m_oldSize - m_oldNumDeleted + 1) / 2;
        rehash_status = REHASH_STATUS::THREE_QUARTER; //updating rehash status
    }
    else {
        datapoints = m_oldSize - m_oldNumDeleted;
        rehash_status = REHASH_STATUS::NOT_REHASHING; //updating rehash status
    }
    for (int i = 0, j = 0; i < datapoints; j++) {
        if (!m_oldTable[j].m_sequence.empty() && m_oldTable[j].m_sequence != DELETEDKEY) {
            unsigned int index = get_index_cur(m_oldTable[j], false);
            m_currentTable[index] = m_oldTable[j];
            m_oldTable[j] = DELETED;
            i++; //only increasing i on transfer
        }
    }
    //size adjustments
    m_currentSize += datapoints;
    m_oldNumDeleted += datapoints;
    if (rehash_status == REHASH_STATUS::NOT_REHASHING) {
        delete[] m_oldTable;
        m_oldTable = nullptr;
        m_oldCap = 0;
        m_oldNumDeleted = 0;
        m_oldSize = 0;
    }
}
