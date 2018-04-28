#include "readpara.h"


ReaderFileAss::ReaderFileAss (string fn1, string fn2, int maxLine) {
   this->scoreFileName = fn1;
   this->targetFileName = fn2;
   this->maxLine = maxLine;
   this->nusr = 0;
   this->nmsg = 0;
}

void ReaderFileAss::getDataDimension () {
   ifstream inf;
   inf.open(this->scoreFileName.c_str(), ifstream::in);
   cout << " read " << this->scoreFileName.c_str() << endl;

   size_t index_separator;
   size_t index_separator_1;
   int iusr, imsg, usrNum, msgNum;
   int k = 0;
   double point;


   // read all usrId and msgId and save them in a set.
   set<string> usrIds;
   set<string> msgIds;
   string usrId, msgId;
   string str, line;


   while (!inf.eof()) {
      k += 1;
      if (k > maxLine) break;

      getline(inf, line);
      //if (k % 1000 == 0)
      //cout << line.c_str() << endl;
      //  cout << " k = " << k << endl;

      index_separator = line.find(',', 0);
      usrId = line.substr(0, index_separator);

      index_separator_1 = line.find(',', index_separator + 1);
      msgId = line.substr(index_separator + 1, index_separator_1 - index_separator - 1);

      if (!usrId.empty()) usrIds.insert(usrId);
      if (!msgId.empty()) msgIds.insert(msgId);

      // cout << usrId.c_str() << " " << msgId.c_str() << endl;
   }
   inf.close();

   usrNum = usrIds.size();
   msgNum = msgIds.size();

   cout << " number of users: " << usrNum << endl;
   cout << " number of messages: " << msgNum << endl;

   k = 0;
   for (set<string>::iterator it = usrIds.begin(); it != usrIds.end(); it++) {
       this->id2usr.insert(pair<int,string>(k, *it));
       this->usr2id.insert(pair<string,int>(*it, k));
       k += 1;
   }
   k = 0;
   for (set<string>::iterator it = msgIds.begin(); it != msgIds.end(); it++) {
       this->id2msg.insert(pair<int,string>(k, *it));
       this->msg2id.insert(pair<string,int>(*it, k));
       k += 1;
   }

   this->nusr = usrNum;
   this->nmsg = msgNum;
   this->ns = usrNum * msgNum;
   this->nm = usrNum * (msgNum + 1);
   this->score = new double[this->ns];
   this->upper = new double[this->nmsg];
   this->bottom = new double[this->nmsg];
   this->target = new double[this->nmsg];
}



void ReaderFileAss::readScoreFromFile () {

   size_t index_separator;
   size_t index_separator_1;
   size_t index_separator_2;
   int iusr, imsg, usrNum, msgNum;
   int k = 0;

   double point;
   string usrId, msgId;
   string str, line;

   ifstream inf;
   inf.open(this->scoreFileName.c_str(), ifstream::in);

   memset(this->score, 0, this->ns * sizeof(double));

   while (!inf.eof()) {
      k += 1;
      if (k > maxLine) break;

      getline(inf, line);

      index_separator = line.find(',', 0);
      usrId = line.substr(0, index_separator);

      index_separator_1 = line.find(',', index_separator + 1);
      msgId = line.substr(index_separator + 1, index_separator_1 - index_separator - 1);

      index_separator_2 = line.find(',', index_separator_1 + 1);
      str = line.substr(index_separator_1 + 1, index_separator_2 - index_separator_1 - 1);

      if (!usrId.empty() && !msgId.empty()) {
	  iusr = usr2id[usrId];
          imsg = msg2id[msgId];
          int ks = iusr * nmsg + imsg;
	  if (ks < 0 || ks >= nusr * nmsg) continue;
          score[ks] = atof(str.c_str());
	  //cout << " ks = " << ks << ", score = " << score[ks] << endl;
      }

      // cout << usrId.c_str() << " " << msgId.c_str() << endl;
   }
   inf.close();


   // normalize
   double sm = 0.0;
   for (int i = 0; i < this->ns; i ++) {
       if (sm < score[i]) sm = score[i];
   }
   for (int i = 0; i < this->ns; i ++) {
       score[i] /= sm;
   }

}

void ReaderFileAss::readTargetFromFile () {

   size_t index_sep;
   size_t index_sep1;
   size_t index_sep2;
   size_t index_sep3;
   int iusr, imsg, usrNum, msgNum;
   int k = 0;
   string str, line, msgId;

   ifstream inf;
   inf.open(this->targetFileName.c_str(), ifstream::in);

   for (int imsg = 0; imsg < nmsg; imsg ++) {
       this->target[imsg] = nusr / nmsg;
       this->bottom[imsg] = 0.0;
       this->upper[imsg] = nusr;
   }

   while (!inf.eof()) {
      k += 1;
      if (k > nmsg) break;

      getline(inf, line);

      index_sep = line.find(',', 0);
      msgId = line.substr(0, index_sep);

      if (msgId.empty()) continue;
      imsg = msg2id[msgId];

      index_sep1 = line.find(',', index_sep + 1);
      str = line.substr(index_sep + 1, index_sep1 - index_sep - 1);
      this->bottom[imsg] = atof(str.c_str());

      index_sep2 = line.find(',', index_sep1 + 1);
      str = line.substr(index_sep1 + 1, index_sep2 - index_sep1 - 1);
      this->upper[imsg] = atof(str.c_str());

      index_sep3 = line.find(',', index_sep2 + 1);
      str = line.substr(index_sep2 + 1, index_sep3 - index_sep2 - 1);
      this->target[imsg] = atof(str.c_str());

   }
   inf.close();
}

void ReaderFileAss::setTarget (int flag) {
   if (flag == 0) {
      for (int imsg = 0; imsg < nmsg; imsg ++) {
          this->target[imsg] = nusr / nmsg;
          this->bottom[imsg] = 300.0;
          this->upper[imsg] = nusr;
      }
   }
   else {
      srand((unsigned)time(NULL));
      for (int i = 0; i < nmsg; i ++)
          this->bottom[i] = 200.0 + 300.0 * rand()/(RAND_MAX + 1.0);
   }
}

void ReaderFileAss::printScore () {
   for (int iusr = 0; iusr < nusr; iusr ++)
   for (int imsg = 0; imsg < nmsg; imsg ++) {
       cout << " score[" << iusr << "][" << imsg << "] = " << score[iusr * nmsg + imsg] << endl;
   }
}


void ReaderFileAss::printTarget () {
   double msg_score;
   cout << " imsg " << " msgId " << " score " << "  bottom " << " upper " << " target " << endl;
   for (int imsg = 0; imsg < nmsg; imsg ++) {
       msg_score = 0;	
       for (int iusr = 0; iusr < nusr; iusr ++) 
	   msg_score += score[iusr * nmsg + imsg];
       cout << imsg << " " << id2msg[imsg].c_str() << " " << msg_score << " " << bottom[imsg] << " " << upper[imsg] << " " << target[imsg] << endl;
   }
}

void ReaderFileAss::free() {
   delete [] this->score;
   delete [] this->upper;
   delete [] this->bottom;
   delete [] this->target;
}
