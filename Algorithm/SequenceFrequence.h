#ifndef SEQUENCEFREQUENCE_H
#define SEQUENCEFREQUENCE_H

#define MIN_FREQUENCE 6
#define FWD 0
#define RC 1
#define MAX_FREQUENCE 200


class SequenceFrequence
{
	public:
		SequenceFrequence(std::string, ErrorCorrectParameters& m_params);
		~SequenceFrequence();
		size_t isFrequenceEnough(ErrorCorrectParameters& m_params);
		bool isCompletelyCorrect(std::string &read_sequence, ErrorCorrectParameters& m_params);
		void outputKmerFrequence(ErrorCorrectParameters& m_params);
		void refreshFrequence(std::string read_sequence, ErrorCorrectParameters& m_params);
		int isFrequenceChange(std::string &read_sequence, ErrorCorrectParameters& m_params);
		//void setIdentityVector(int round, int identity);
		//void outputIdentity(int round);
		void outputFrequnecInterval(int round);
		bool isKmerNeedDecrease(ErrorCorrectParameters& m_params);
		void isKmerNeedIncrease(int round);
		void increaseKmer();
		void decreaseKmer();
		size_t continueInterval();
		size_t getklen();
		std::string getquery();
		int isKmerchange();
		size_t getkmerLength();
		
	private:
		std::string query;
		std::vector<std::vector <long> > queryFrequence;
		//std::vector<std::vector <long> > identityVector;
		std::vector<std::vector <long> > intervalVector;
		size_t invervalIndex;
		size_t klen;
		size_t kmerLength;
		size_t changeKmer;
};

#endif