#include "SequenceFrequence.h"
#include "ErrorCorrectProcess.h"
#include <math.h>
#include <vector>
#include <iostream>

SequenceFrequence::SequenceFrequence(std::string current_sequence, ErrorCorrectParameters& m_params)
{
	size_t fwd_count, rc_count;
	query = current_sequence;
	klen = current_sequence.size() - m_params.kmerLength + 1;
	kmerLength = m_params.kmerLength;
	changeKmer = 0;
	invervalIndex = 0;
	
	queryFrequence.resize(2);
	intervalVector.resize(m_params.numOverlapRounds);
	//identityVector.resize(m_params.numOverlapRounds);
	for(int i = 0; i < 2; i++)
		queryFrequence[i].resize(klen);
	for(int i = 0; i < klen; i++)
	{
			//std::cout<<"i:\t"<<i<<std::endl;
			std::string kmer = query.substr(i, kmerLength);
			fwd_count = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
			rc_count = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer),m_params.indices);
			queryFrequence[FWD][i] = fwd_count;
			queryFrequence[RC][i] = rc_count;
	}
	for(int i = 0; i < m_params.numOverlapRounds; i++)
	{
		//identityVector[i].resize(100);
		intervalVector[i].resize(klen);
	}
}

//
SequenceFrequence::~SequenceFrequence()
{

}

//
size_t SequenceFrequence::isFrequenceEnough(ErrorCorrectParameters& m_params)
{
	size_t i, fwd_count,rc_count, flag = 1,count = 0;
	//int nk = readSequence.size()-m_params.kmerLength;
	//	std::cout<<"readSequence.size=="<<readSequence.size()<<" m_params.kmerLength=="<<m_params.kmerLength<<"\n";
	for(i = 0; i < klen; i++)
	{
		std::string kmer = query.substr(i, kmerLength);
		fwd_count = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
		rc_count = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer),m_params.indices);
		queryFrequence[FWD][i] = fwd_count;
		queryFrequence[RC][i] = rc_count;
		if(fwd_count > MIN_FREQUENCE && rc_count > MIN_FREQUENCE)
			continue;
		else
		{
			count++;
			flag = 0;
		}
	}
	std::cout<<"count == "<<count<<"klen == "<<klen<<"\n";
	//std::cout<<"k == "<<m_params.kmerLength<<"\n";
	return flag;
}

//
void SequenceFrequence::refreshFrequence(std::string read_sequence ,ErrorCorrectParameters& m_params)
{
	size_t i, fwd_count,rc_count;
	query.erase();
	query = read_sequence;
	klen = read_sequence.size() - kmerLength + 1;
	for(i = 0; i < 2; i++)
	{
		queryFrequence[i].clear();
		queryFrequence[i].resize(klen);
	}
	for(i = 0; i < klen; i++)
	{
		std::string kmer = query.substr(i, kmerLength);
		fwd_count = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
		rc_count = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer),m_params.indices);
		queryFrequence[FWD][i] = fwd_count;
		queryFrequence[RC][i] = rc_count;
		//std::cout<<"klen == ";
		//std::cout<<klen<<":"<<fwd_count<<":"<<rc_count<<"\n";
		//std::cout<<i<<":"<<queryFrequence[FWD][i]<<":"<<queryFrequence[RC][i]<<"\n";
	}
}

//
void SequenceFrequence::outputKmerFrequence(ErrorCorrectParameters& m_params)
{
	size_t i, fwd_count,rc_count;
	//int nk = readSequence.size()-m_params.kmerLength;
	//std::cout<<"kmerLength == "<<kmerLength<<"\n";
	for(i = 0; i < klen; i++)
	{
		std::string kmer = query.substr(i, kmerLength);
		fwd_count = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
		rc_count = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer),m_params.indices);
		std::cout<<i<<":"<<fwd_count<<":"<<rc_count<<"\n";

	}
}

//
bool SequenceFrequence::isCompletelyCorrect(std::string &read_sequence, ErrorCorrectParameters& m_params)
{
	bool flag = true;
	size_t i, fwd_count, rc_count;
	if(queryFrequence[FWD].size() == read_sequence.size())
	{
		for(i = 0; i < klen; i++)
		{
			std::string kmer = query.substr(i, kmerLength);
			fwd_count = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
			rc_count = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer),m_params.indices);
			if(fwd_count <= MIN_FREQUENCE && rc_count <= MIN_FREQUENCE)
				flag=false;
			if(queryFrequence[FWD][i] != fwd_count || queryFrequence[RC][i] != rc_count)
			{
				queryFrequence[FWD][i] = fwd_count;
				queryFrequence[RC][i] = rc_count;
				//std::cout<<"frequence not the same\n";
			}
		}
	}	
	else
	{   
		query.erase();
		query = read_sequence;
		klen = read_sequence.size() - kmerLength + 1;
		
		for(i = 0; i < 2; i++)
		{
			queryFrequence[i].clear();
			queryFrequence[i].resize(klen);
		}
		
		for(i = 0; i < klen; i++)
		{
			//std::cout<<"i:\t"<<i<<std::endl;
			std::string kmer = query.substr(i, kmerLength);
			fwd_count = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
			rc_count = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer),m_params.indices);
			queryFrequence[FWD][i] = fwd_count;
			queryFrequence[RC][i] = rc_count;
			if(fwd_count < MIN_FREQUENCE && fwd_count < MIN_FREQUENCE)
			{
				flag = false;
			}
		}
	}
	return flag;
}

//
int SequenceFrequence::isFrequenceChange(std::string &read_sequence, ErrorCorrectParameters& m_params)
{
	size_t i,fwd_count = 0, rc_count = 0;
	for(i = 0; i < klen; i++)
	{
		std::string kmer = read_sequence.substr(i, m_params.kmerLength);
		fwd_count = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
		rc_count = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer),m_params.indices);

		if(queryFrequence[FWD][i] != fwd_count || queryFrequence[RC][i] != rc_count)
			return 1;
	}

	return 0;
}

/*void SequenceFrequence::setIdentityVector(int round, int identity)
{
	identityVector[round][identity]++;
}
void SequenceFrequence::outputIdentity(int round)
{
	double mean = 0.0, temp_mean = 0.0,temp = 0.0;
	for(int i = 0; i < 100; i++)
	{	//count*identity
		mean+=identityVector[round][i]*i;
		temp+=identityVector[round][i];
	}
	mean=mean/temp;
	
	for(int i = 0; i < 100; i++)
		//count*identity^2
		temp_mean+=identityVector[round][i]*pow(i,2);
	std::cout<<"-----------outputIdentity------------"<<std::endl;
	std::cout<<"\tround "<<round;
	std::cout<<"\tmean identity :\t"<<mean<<std::endl;
	std::cout<<"\tSD identity :\t"<<sqrt(temp_mean/temp - pow(mean,2))<<std::endl;
	std::cout<<"-------------------------------------"<<std::endl;
}
*/

//
void SequenceFrequence::outputFrequnecInterval(int round)
{
	std::cout<<"-------outputFrequnecInterval--------"<<std::endl;
	std::cout<<"\tround :\t"<<round<<std::endl;
	for(size_t i = 0; i < invervalIndex;i++)
		std::cout<<"\tfrequence interval :\t"<<intervalVector[round][i]<<std::endl;
	invervalIndex = 0;
	std::cout <<"-------------------------------------"<<std::endl;
}

//
bool SequenceFrequence::isKmerNeedDecrease(ErrorCorrectParameters& m_params)
{
	size_t i,count = 0;
	long fwd_count,rc_count;
	std::cout<<"--------in isKmerNeedDecrease--------\n";
	std::cout<<"\toriginal kmer :\t"<<kmerLength<<"\n";
	if(kmerLength - 3 <= 0)
		return false;
	else
	{
		for(i = 0; i < klen; i++)
			if(queryFrequence[FWD][i] <= 10 && queryFrequence[RC][i] <= 10)
				count++;
		
		//std::cout<<"count == "<<count<<"klen == "<<klen<<"\n";
		//std::cout<<"kmerLength == "<<kmerLength<<"\n";
	}

	if(count == klen)
		return true;
	else 
		return false;
}

//
void SequenceFrequence::isKmerNeedIncrease(int round)
{
	size_t i, temp_count, max_count = 0, total_interval = 0;
	invervalIndex = 0;
	std::cout<<"--------in isKmerNeedIncrease--------\n";
	if(kmerLength+3>=100)
		kmerLength = kmerLength;
	else
	
		for(i = 0; i < klen; i++)
		{
			if(queryFrequence[FWD][i] > MIN_FREQUENCE || queryFrequence[RC][i] > MIN_FREQUENCE)
			{
				temp_count = 0;
				for(;i < klen; i++)
				{
					if(queryFrequence[FWD][i] > MIN_FREQUENCE || queryFrequence[RC][i] > MIN_FREQUENCE)
						temp_count++;
					else
						break;
				}
				std::cout<<"\tcontinue inverval :\t" <<temp_count<<"\n";
				intervalVector[round][invervalIndex] = temp_count;
				invervalIndex++;
				total_interval += temp_count + kmerLength - 1;
				if(temp_count > max_count && temp_count < 100 - kmerLength +1)
					max_count = temp_count;
			}
		}
		
	std::cout<<"\tquery length :\t"<<query.size()<<"\n";
	std::cout<<"\ttotal_interval :\t"<<total_interval<<"\n";
	std::cout<<"\toriginal kmer :\t"<<kmerLength<<"\n";
	std::cout<<"\tmax_count :\t"<<max_count<<"\n";
	
	if(kmerLength + max_count - 1 > kmerLength)
		kmerLength += max_count - 1;
	
	klen = query.size() - kmerLength + 1;
	std::cout<<"\tincrease kmer :\t"<<kmerLength<<"\n";
	std::cout <<"-------------------------------------"<<std::endl;
}

//
void SequenceFrequence::decreaseKmer()
{
	kmerLength -= 2;
	klen = query.size() - kmerLength + 1;
	changeKmer++;
	std::cout<<"\tdecrease kmer :\t"<<kmerLength<<"\n";
	std::cout <<"-------------------------------------"<<std::endl;
}

//
void SequenceFrequence::increaseKmer()
{
	kmerLength += 3;
	klen = query.size() - kmerLength + 1;
	changeKmer--;
	std::cout<<"increase kmer to "<<kmerLength<<"\n";
}

size_t SequenceFrequence::continueInterval()
{
	size_t count = 0, sumInterval = 0, sumCount = 0, sumSeedNum = 0;
	size_t countTimes = 0, countMax = 0,countMin = 100000;
	std::cout<<"\tkmerLength : "<<kmerLength<<std::endl;
	for(size_t i = 0; i < klen; i++)
	{
		if((queryFrequence[FWD][i] >= 6 || queryFrequence[RC][i] >= 6) && count == 0)
		{
			for(; i < klen; i++)
			{
				if(queryFrequence[FWD][i] >= 6 || queryFrequence[RC][i] >= 6)
					count++;
				else 
					break;
			}
		}
		else if(count != 0)
		{
			countTimes++;
			sumCount += count;
			sumInterval += count - 1 + kmerLength;
			sumSeedNum += count*20*2;
			//std::cout<<"\tcontinueInterval : "<<count - 1 + kmerLength<<std::endl;
			if(countMax < count - 1 + kmerLength)
				countMax = count - 1 + kmerLength;
			if(countMin > count - 1 + kmerLength)
				countMin = count - 1 + kmerLength;
			count = 0;
		}
		
	}

	std::cout<<"\tcountTimes : "<<countTimes<<std::endl;
	std::cout<<"\tcountMin : "<<countMin<<std::endl;
	std::cout<<"\tcountAve : "<<(double)sumInterval/countTimes<<std::endl;
	std::cout<<"\tcountMax : "<<countMax<<std::endl;
	std::cout<<"\tsumCount : "<<sumCount<<std::endl;
	std::cout<<"\tsumSeedNum : "<<sumSeedNum<<std::endl;
	
	return sumInterval;
}

//
size_t SequenceFrequence::getklen()
{
	return klen;
}

//
std::string SequenceFrequence::getquery()
{
	return query;
}

//
int SequenceFrequence::isKmerchange()
{
	return changeKmer;
}

//
size_t SequenceFrequence::getkmerLength()
{
	return kmerLength;
}