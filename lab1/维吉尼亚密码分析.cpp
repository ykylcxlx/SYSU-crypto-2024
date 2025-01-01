#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <unordered_map>
#include <algorithm>
#include <numeric>
using namespace std;

double getSum(const vector<double>& a){
    double sum = 0;
    for(const auto& val : a){
        sum += val;
    }
    return sum;
}

int getKeyLen(const string& low_text){
    int temp_len = low_text.size() / 50;
    vector<int> v(temp_len - 1);
    for(int i = 2; i <= temp_len; i++){
        v[i - 2] = i;
    }
    double max_mean = 0;
    int key_len = 0;
    for(const int& num : v){
        vector<vector<double>> freq(num, vector<double>(26, 0));

        for(int j = 0; j < low_text.size(); j += 1){
            freq[j % num][low_text[j] - 'a'] += 1;
        }

        for(auto& row : freq){
            double sum = getSum(row);
            for(auto& val : row){
                val /= sum;
                val *= val;
            }
        }

        double mean = 0;
        for(const auto& row : freq){
            mean += getSum(row);
        }
        mean /= freq.size();

        if(max_mean < mean){
            max_mean = mean;
            key_len = num;
        }
    }
    return key_len;
}

void computeLPSArray(const string& pat, int M, vector<int>& lps) {
    int len = 0;
    lps[0] = 0;

    int i = 1;
    while (i < M) {
        if (pat[i] == pat[len]) {
            len++;
            lps[i] = len;
            i++;
        } else {
            if (len != 0) {
                len = lps[len - 1];
            } else {
                lps[i] = 0;
                i++;
            }
        }
    }
}

string findSmallestRepeatingSubstring(const string& str) {
    if(str.empty()){
        return str;
    }
    int N = str.length();
    vector<int> lps(N);
    computeLPSArray(str, N, lps);
    int len = lps[N - 1];

    if (len == 0 || N % (N - len) != 0) {
        return str;
    }

    return str.substr(0, N - len);
}

int main(){
    string text, line;
    double pr[26] = {0.0726, 0.0027, 0.0403, 0.0403, 0.1371, 0.0188, 0.0323, 0.0349, 0.0806, 0.0054, 0.0027, 0.0403, 0.0403, 0.0968, 0.0565, 0.0215, 0.0027, 0.0565, 0.0538, 0.0968, 0.0403, 0.0054, 0.0027, 0.0027, 0.0108, 0.0054};

    while(cin >> line){
        text += line + ' ';
    }

    string low_text = text;
    low_text.erase(remove_if(low_text.begin(), low_text.end(), [](char c) { return !isalpha(c); }), low_text.end());
    
    transform(low_text.begin(), low_text.end(), low_text.begin(), ::tolower);

    int key_len = getKeyLen(low_text);

    string key;
    for(int i = 0; i < key_len; i++){
        double max_Mg = 0;
        int id = 0;
        vector<double> f(26, 0);

        for(int j = i; j < low_text.size(); j += key_len){
            f[low_text[j] - 'a'] += 1;
        }

        double sumj = getSum(f);

        for(auto& val : f){
            val /= sumj;
        }

        for(int g = 0; g < 26; g++){
            double Mg = 0;
            for(int i = 0; i < 26; i++){
                Mg += pr[i] * f[(i + g) % 26];
            }
            if(Mg > max_Mg){
                max_Mg = Mg;
                id = g;
            }
        }
        key.push_back('A' + id);
    }

    string keyy = findSmallestRepeatingSubstring(key);

    for(const char& c : keyy){
        cout << c;
    }
    cout << endl;

    if(key.empty()){
        cout << "can't find key" << endl;
        return 0;
    }

    int j = 0;
    for(int i = 0,j=0;i<text.size();i++){
        
        if(text[i]>='a'&&text[i]<='z')
        text[i] = ((low_text[j]+26-(key[j%key_len]-'A'))%97)%26+97;
        else if(text[i]>='A'&&text[i]<='Z')
        text[i] = ((low_text[j]+26-(key[j%key_len]-'A'))%97)%26+65;
        else
        continue;

        j++;
        
    }

    cout << text;
    return 0;
}