#include "dnadb.h"
#include <random>
#include <vector>
#include <set>
#include <algorithm>
enum RANDOM { UNIFORMINT, UNIFORMREAL, NORMAL };
class Random {
public:
    Random(int min, int max, RANDOM type = UNIFORMINT, int mean = 50, int stdev = 20) : m_min(min), m_max(max), m_type(type)
    {
        if (type == NORMAL) {
            //the case of NORMAL to generate integer numbers with normal distribution
            m_generator = std::mt19937(m_device());
            //the data set will have the mean of 50 (default) and standard deviation of 20 (default)
            //the mean and standard deviation can change by passing new values to constructor
            m_normdist = std::normal_distribution<>(mean, stdev);
        }
        else if (type == UNIFORMINT) {
            //the case of UNIFORMINT to generate integer numbers
            // Using a fixed seed value generates always the same sequence
            // of pseudorandom numbers, e.g. reproducing scientific experiments
            // here it helps us with testing since the same sequence repeats
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_unidist = std::uniform_int_distribution<>(min, max);
        }
        else { //the case of UNIFORMREAL to generate real numbers
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_uniReal = std::uniform_real_distribution<double>((double)min, (double)max);
        }
    }
    void setSeed(int seedNum) {
        // we have set a default value for seed in constructor
        // we can change the seed by calling this function after constructor call
        // this gives us more randomness
        m_generator = std::mt19937(seedNum);
    }

    int getRandNum() {
        // this function returns integer numbers
        // the object must have been initialized to generate integers
        int result = 0;
        if (m_type == NORMAL) {
            //returns a random number in a set with normal distribution
            //we limit random numbers by the min and max values
            result = m_min - 1;
            while (result < m_min || result > m_max)
                result = m_normdist(m_generator);
        }
        else if (m_type == UNIFORMINT) {
            //this will generate a random number between min and max values
            result = m_unidist(m_generator);
        }
        return result;
    }

    double getRealRandNum() {
        // this function returns real numbers
        // the object must have been initialized to generate real numbers
        double result = m_uniReal(m_generator);
        // a trick to return numbers only with two deciaml points
        // for example if result is 15.0378, function returns 15.03
        // to round up we can use ceil function instead of floor
        result = std::floor(result * 100.0) / 100.0;
        return result;
    }

private:
    int m_min;
    int m_max;
    RANDOM m_type;
    std::random_device m_device;
    std::mt19937 m_generator;
    std::normal_distribution<> m_normdist;//normal distribution
    std::uniform_int_distribution<> m_unidist;//integer uniform distribution
    std::uniform_real_distribution<double> m_uniReal;//real uniform distribution

};
class Tester {
public:
    bool test_insert();
    bool test_find();
    bool test_find_colliding();
    bool test_remove();
    bool test_remove_colliding();
    bool test_rehash_insertion();
    bool test_rehash_removal();
};

unsigned int hashCode(const string str);
string sequencer(int size, int seedNum);

int main() {
    Tester tester;
    tester.test_insert();
    cout << endl;
    tester.test_find();
    cout << endl;
    tester.test_find_colliding();
    cout << endl;
    tester.test_remove();
    cout << endl;
    tester.test_remove_colliding();
    cout << endl;
    tester.test_rehash_insertion();
    cout << endl;
    tester.test_rehash_removal();
    return 0;
}
unsigned int hashCode(const string str) {
    unsigned int val = 0;
    const unsigned int thirtyThree = 33;  // magic number from textbook
    for (int i = 0; i < str.length(); i++)
        val = val * thirtyThree + str[i];
    return val;
}
string sequencer(int size, int seedNum) {
    //this function returns a random DNA sequence
    string sequence = "";
    Random rndObject(0, 3);
    rndObject.setSeed(seedNum);
    for (int i = 0; i < size; i++) {
        sequence = sequence + ALPHA[rndObject.getRandNum()];
    }
    return sequence;
}

bool Tester::test_insert() {
    cout << "Testing Insert Function:" << endl;
    set<unsigned int> used_hash_indices;
    vector<DNA> dataList;
    Random RndLocation(MINLOCID, MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);
    bool result = true;
    for (int i = 0; i < 49; i++) {
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        // saving data for later use
        unsigned int idx = hashCode(dataObj.getSequence()) % dnadb.m_currentCap;
        if (used_hash_indices.count(idx) == 0) {
            used_hash_indices.insert(idx);
            dataList.push_back(dataObj);// Only retain the DNAs that hash to different
        }                                   // location (to ensure no collisions when insert-
        // ing into the hash table).
    }
    // inserting data in to the DnaDb object
    cout << "Inserting " << dataList.size() << " non-colliding DNA Objects." << endl;
    for (const auto& D : dataList) {
        dnadb.insert(D);    // Shouldn't be any collisions here

        // Check if the item was inserted at the correct index.
        unsigned int idx = hashCode(D.getSequence()) % dnadb.m_currentCap;
        if (!(dnadb.m_currentTable[idx] == D)) {
            cout << "Insert Failed!" << endl;
            return false;
        }
    }

    // Check that the correct number of elements was inserted
    if (dataList.size() != dnadb.m_currentSize) {
        return false;
    }
    cout << "Insert Succeeded!" << endl;
    return true;
}

bool Tester::test_find() {
    cout << "Testing Find Error with Empty Table: " << endl;
    DnaDb dnadb(MINPRIME, hashCode);
    DNA test = dnadb.getDNA("AAGCT", 21);
    if (test == EMPTY) {
        cout << "Test Successful!" << endl;
    }
    else {
        cout << "Test Failed" << endl;
        return false;
    }
    cout << "Testing Find Error with Full Table: " << endl;
    set<unsigned int> used_hash_indices;
    vector<DNA> dataList;
    Random RndLocation(MINLOCID, MAXLOCID);
    for (int i = 0; i < 49; i++) {
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        // saving data for later use
        unsigned int idx = hashCode(dataObj.getSequence()) % dnadb.m_currentCap;
        if (used_hash_indices.count(idx) == 0) {
            used_hash_indices.insert(idx);
            dataList.push_back(dataObj);// Only retain the DNAs that hash to different
        }
    }
    cout << "Inserting " << dataList.size() << " non-colliding DNA Objects." << endl;
    for (const auto& D : dataList) {
        dnadb.insert(D);    // Shouldn't be any collisions here
    }
    int i = 0;
    DNA new_DNA;
    while (true) {
        new_DNA = DNA(sequencer(5, i++), RndLocation.getRandNum());
        if (std::find(dataList.cbegin(), dataList.cend(), new_DNA) == dataList.cend()) {
            break;
        }
    }
    test = dnadb.getDNA(new_DNA.getSequence(), new_DNA.getLocId());
    if (test == EMPTY) {
        cout << "Test Successful!" << endl;
    }
    else {
        cout << "Test Failed" << endl;
        return false;
    }
    new_DNA = dataList.front();
    cout << "Testing Find with a full table: " << endl;
    test = dnadb.getDNA(new_DNA.getSequence(), new_DNA.getLocId());
    if (test == EMPTY) {
        cout << "Test Failed!" << endl;
        return false;
    }
    else {
        cout << "Test Successful" << endl;
    }
    new_DNA = dataList.back();
    cout << "Testing Find with a full table: " << endl;
    test = dnadb.getDNA(new_DNA.getSequence(), new_DNA.getLocId());
    if (test == EMPTY) {
        cout << "Test Failed!" << endl;
        return false;
    }
    else {
        cout << "Test Successful" << endl;
    }
    return true;
}

bool Tester::test_find_colliding() {
    cout << "Testing Find with Colliding Data" << endl;
    DnaDb dnadb(MINPRIME, hashCode);
    vector<DNA> dataList;
    Random RndLocation(MINLOCID, MAXLOCID);
    for (int i = 0; i < 49; i++) {
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        if (std::find(dataList.cbegin(), dataList.cend(), dataObj) == dataList.cend()) {
            dataList.push_back(dataObj);// Only retain the DNAs that hash to different
        }
    }
    cout << "Inserting " << dataList.size() << " colliding DNA Objects." << endl;
    for (const auto& D : dataList) {
        dnadb.insert(D);
    }

    //Looping through vector
    for (const auto& D : dataList) {
        DNA temp = dnadb.getDNA(D.getSequence(), D.getLocId());
        if (temp == EMPTY) {
            cout << "Test Failed!" << endl;
            return false;
        }
    }
    cout << "Test Successful" << endl;
    return true;
}

bool Tester::test_remove(){
    cout << "Testing Remove Function:" << endl;
    set<unsigned int> used_hash_indices;
    vector<DNA> dataList;
    Random RndLocation(MINLOCID, MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);
    bool result = true;
    for (int i = 0; i < 49; i++) {
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        // saving data for later use
        unsigned int idx = hashCode(dataObj.getSequence()) % dnadb.m_currentCap;
        if (used_hash_indices.count(idx) == 0) {
            used_hash_indices.insert(idx);
            dataList.push_back(dataObj);// Only retain the DNAs that hash to different
        }                                   // location (to ensure no collisions when insert-
        // ing into the hash table).
    }
    // inserting data in to the DnaDb object
    cout << "Inserting " << dataList.size() << " non-colliding DNA Objects." << endl;
    for (const auto& D : dataList) {
        dnadb.insert(D);    // Shouldn't be any collisions here
    }
    cout << "Removing all Objects" << endl;
    for (const auto& D : dataList) {
        if (dnadb.remove(D) == false){
            cout << "Operation Failed" << endl;
            return false;
        }
    }
    if (dnadb.m_currentSize != 0){
        cout << "Operation Failed" << endl;
        return false;
    }
    cout << "Operation Successful" << endl;
    return true;
}

bool Tester::test_remove_colliding(){
    cout << "Testing Remove with Colliding Data" << endl;
    DnaDb dnadb(MINPRIME, hashCode);
    vector<DNA> dataList;
    Random RndLocation(MINLOCID, MAXLOCID);
    for (int i = 0; i < 49; i++) {
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        if (std::find(dataList.cbegin(), dataList.cend(), dataObj) == dataList.cend()) {
            dataList.push_back(dataObj);// Only retain the DNAs that hash to different
        }
    }
    cout << "Inserting " << dataList.size() << " colliding DNA Objects." << endl;
    for (const auto& D : dataList) {
        dnadb.insert(D);
    }
    cout << "Removing all Objects" << endl;
    for (const auto& D : dataList) {
        if (dnadb.remove(D) == false){
            cout << "Operation Failed" << endl;
            return false;
        }
    }
    if (dnadb.m_currentSize != 0){
        cout << "Operation Failed" << endl;
        return false;
    }
    cout << "Operation Successful" << endl;
    return true;
}

bool Tester::test_rehash_insertion(){
    cout << endl << "testing Rehashing During Insertions" << endl;
    DnaDb dnadb(MINPRIME, hashCode);
    vector<DNA> dataList;
    Random RndLocation(MINLOCID, MAXLOCID);
    for (int i = 0; i < 99; i++){
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        if (std::find(dataList.cbegin(), dataList.cend(), dataObj) == dataList.cend()){
            dataList.push_back(dataObj);// Only retain the DNAs that hash to different
        }
    }
    cout << "Current Table Size: " << dnadb.m_currentSize << endl;
    cout << "Current Table Capacity: " << dnadb.m_currentCap << endl;
    cout << "Old Table Size: " << dnadb.m_oldSize << endl;
    cout << "Old Table Capacity: " << dnadb.m_oldCap << endl;
    cout << "Inserting " << dataList.size() << " colliding DNA Objects." << endl << endl;
    for (const auto& D : dataList){
        dnadb.insert(D);
        if (dnadb.rehash_status != DnaDb::REHASH_STATUS::NOT_REHASHING){
            if (dnadb.rehash_status == DnaDb::REHASH_STATUS::QUARTER){
                cout << "Rehash Triggered!" << endl;
            }
            else{
                cout << "Rehash still in Progress!" << endl;
            }
            cout << "Current Table Size: " << dnadb.m_currentSize << endl;
            cout << "Current Table Capacity: " << dnadb.m_currentCap << endl;
            cout << "Old Table Size: " << dnadb.m_oldSize << endl;
            cout << "Old Table Deleted Size: " << dnadb.m_oldNumDeleted << endl;
            cout << "Old Table Capacity: " << dnadb.m_oldCap << endl << endl;
        }
    }
    cout << endl << "Insertions Completed!" << endl;
    cout << "Current Table Size: " << dnadb.m_currentSize << endl;
    cout << "Current Table Capacity: " << dnadb.m_currentCap << endl;
    cout << "Old Table Size: " << dnadb.m_oldSize << endl;
    cout << "Old Table Capacity: " << dnadb.m_oldCap << endl;
    return true;
}

bool Tester::test_rehash_removal() {
    cout << endl << "Testing Rehashing During Removals" << endl;
    DnaDb dnadb(MINPRIME, hashCode);
    vector<DNA> dataList;
    Random RndLocation(MINLOCID, MAXLOCID);
    for (int i = 0; i < 99; i++) {
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        if (std::find(dataList.cbegin(), dataList.cend(), dataObj) == dataList.cend()) {
            dataList.push_back(dataObj);// Only retain the DNAs that hash to different
        }
    }
    cout << "Inserting " << dataList.size() << " colliding DNA Objects." << endl;
    for (const auto& D : dataList) {
        dnadb.insert(D);
    }
    cout << endl << "Insertions Completed!" << endl;
    cout << "Current Table Size: " << dnadb.m_currentSize << endl;
    cout << "Current Table Capacity: " << dnadb.m_currentCap << endl;
    cout << "Old Table Size: " << dnadb.m_oldSize << endl;
    cout << "Old Table Capacity: " << dnadb.m_oldCap << endl << endl;
    cout << "Removing all Objects" << endl;
    for (int i = 0, I = dataList.size(); i < I; i++){ //looping through everything
        DNA D = dataList.at(i);
        if (dnadb.remove(D) == false){
            cout << "Operation Failed" << endl;
            return false;
        }
        if (dnadb.rehash_status != DnaDb::REHASH_STATUS::NOT_REHASHING){
            if (dnadb.rehash_status == DnaDb::REHASH_STATUS::QUARTER){
                cout << "Rehash Triggered!" << endl;
            }
            else{
                cout << "Rehash still in Progress!" << endl;
            }
            cout << "Current Table Size: " << dnadb.m_currentSize << endl;
            cout << "Current Table Deleted Size: " << dnadb.m_currNumDeleted << endl;
            cout << "Current Table Capacity: " << dnadb.m_currentCap << endl;
            cout << "Old Table Size: " << dnadb.m_oldSize << endl;
            cout << "Old Table Deleted Size: " << dnadb.m_oldNumDeleted << endl;
            cout << "Old Table Capacity: " << dnadb.m_oldCap << endl << endl;
        }
    }
    cout << endl << "Removals Completed!" << endl;
    cout << "Current Table Size: " << dnadb.m_currentSize << endl;
    cout << "Current Table Deleted Size: " << dnadb.m_currNumDeleted << endl;
    cout << "Current Table Capacity: " << dnadb.m_currentCap << endl;
    cout << "Old Table Size: " << dnadb.m_oldSize << endl;
    cout << "Old Table Deleted Size: " << dnadb.m_oldNumDeleted << endl;
    cout << "Old Table Capacity: " << dnadb.m_oldCap << endl << endl;
    return true;
}