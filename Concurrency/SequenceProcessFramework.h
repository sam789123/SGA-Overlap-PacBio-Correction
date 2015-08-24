//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SequenceProcessFramework - Generic framework for performing
// some operations on input data produced by a generator,
// serially or in parallel. 
//
#include "ThreadWorker.h"
#include "Timer.h"
#include "SequenceWorkItem.h"
#include "config.h"

#if HAVE_OPENMP
#include <omp.h>
#endif

#ifndef SEQUENCEPROCESSFRAMEWORK_H
#define SEQUENCEPROCESSFRAMEWORK_H

namespace SequenceProcessFramework
{

const size_t BUFFER_SIZE = 100;

// Generic function to process n work items from a file. 
// With the default value of -1, n becomes the largest value representable for
// a size_t and all values will be read
template<class Input, class Output, class Generator, class Processor, class PostProcessor>
size_t processWorkSerial(Generator& generator, Processor* pProcessor, PostProcessor* pPostProcessor, size_t n = -1)
{
    Timer timer("SequenceProcess", true);
    Input workItem;
    
    // Generate work items using the generic generation class while the number
    // of sequences consumed from the SeqReader is less than n and there 
    // are still sequences to consume from the reader
    while(generator.getNumConsumed() < n && generator.generate(workItem))
    {
        Output output = pProcessor->process(workItem);
        
        pPostProcessor->process(workItem, output);
        if(generator.getNumConsumed() % 50000 == 0)
            printf("[sga] Processed %zu sequences (%lfs elapsed)\n", generator.getNumConsumed(), timer.getElapsedWallTime());
    }

    assert(n == (size_t)-1 || generator.getNumConsumed() == n);

    //
    double proc_time_secs = timer.getElapsedWallTime();
    printf("[sga::process] processed %zu sequences in %lfs (%lf sequences/s)\n", 
            generator.getNumConsumed(), proc_time_secs, (double)generator.getNumConsumed() / proc_time_secs);    
    
    return generator.getNumConsumed();
}

// Wrapper function for performing operations over n elements from a SeqReader
// This wrapper exists for interfacing with legacy code
template<class Input, class Output, class Processor, class PostProcessor>
size_t processSequencesSerial(SeqReader& reader, Processor* pProcessor, PostProcessor* pPostProcessor, size_t n = -1)
{
    WorkItemGenerator<Input> generator(&reader);
    return processWorkSerial<Input, 
                             Output, 
                             WorkItemGenerator<Input>, 
                             Processor, 
                             PostProcessor>(generator, pProcessor, pPostProcessor, n);
}

// Wrapper function for performing operations over every sequence read in readsFile
template<class Input, class Output, class Processor, class PostProcessor>
size_t processSequencesSerial(const std::string& readsFile, Processor* pProcessor, PostProcessor* pPostProcessor)
{
    SeqReader reader(readsFile);
    WorkItemGenerator<Input> generator(&reader);
    return processWorkSerial<Input, 
                             Output, 
                             WorkItemGenerator<Input>, 
                             Processor, 
                             PostProcessor>(generator, pProcessor, pPostProcessor);
}


// Design:
// This function is a generic function to read some INPUT from a 
// generic generator object, then perform work on them.
// The actual processing is done by the Processor class 
// that is passed in. The number of threads
// created is determined by the size of the vector of processors - 
// one thread per processor. 
//
// The function buffers batches of input data.
// Once the buffers are full, the reads are dispatched to the thread
// which run the actual processing independently. An optional post processor
// can be specified to process the results that the threads return. If the n
// parameter is used, at most n sequences will be read from the file.
// 
// This version is based on pthreads.
template<class Input, class Output, class Generator, class Processor, class PostProcessor>
size_t processWorkParallelPthread(Generator& generator, 
                                  std::vector<Processor*> processPtrVector, 
                                  PostProcessor* pPostProcessor, 
                                  size_t n = -1)
{
    Timer timer("SequenceProcess", true);

    // Helpful typedefs
    typedef ThreadWorker<Input, Output, Processor> Thread;
    typedef std::vector<Thread*> ThreadPtrVector;

    typedef std::vector<Input> InputItemVector;
    typedef std::vector<InputItemVector*> InputBufferVector;

    typedef std::vector<Output> OutputVector;
    typedef std::vector<OutputVector*> OutputBufferVector;
    typedef std::vector<sem_t*> SemaphorePtrVector;


    // Initialize threads, one thread per processor that was passed in
    int numThreads = processPtrVector.size();

    ThreadPtrVector threadVec(numThreads);
    InputBufferVector inputBuffers(numThreads);
    OutputBufferVector outputBuffers(numThreads);
    SemaphorePtrVector semVec(numThreads);
	
	
	
    // Create the threads
    for(int i = 0; i < numThreads; ++i)
    {
        semVec[i] = new sem_t;
        int ret = sem_init( semVec[i], PTHREAD_PROCESS_PRIVATE, 0 );
        if(ret != 0)
        {
            std::cerr << "Semaphore initialization failed with error " << ret << "\n";
            std::cerr << "You are probably running on OSX which does not provide unnamed semaphores\n";
            exit(EXIT_FAILURE);
        }

        // Create and start the thread
        threadVec[i] = new Thread(semVec[i], processPtrVector[i], BUFFER_SIZE);
        threadVec[i]->start();

        inputBuffers[i] = new InputItemVector;
        inputBuffers[i]->reserve(BUFFER_SIZE);

        outputBuffers[i] = new OutputVector;
        outputBuffers[i]->reserve(BUFFER_SIZE);
		
    }
	
    size_t numWorkItemsRead = 0;
    size_t numWorkItemsWrote = 0;
    bool done = false;
    int next_thread = 0;
    int num_buffers_full = 0;

    while(!done)
    {
        // Parse reads from the stream and add them into the incoming buffers
        Input workItem;
		
        bool valid = generator.generate(workItem);
		std::cout<<"1---------"<<std::endl;
        if(valid)
        {
            inputBuffers[next_thread]->push_back(workItem);
            numWorkItemsRead += 1;

            // Change buffers if this one is full
            if(inputBuffers[next_thread]->size() == BUFFER_SIZE)
            {
                ++num_buffers_full;
                ++next_thread;
            }
        }
	
        done = !valid || generator.getNumConsumed() == n;
		std::cout<<"2---------"<<std::endl;
        // Once all buffers are full or the input is finished, dispatch the reads to the threads
        // by swapping work buffers. 
        if(num_buffers_full == numThreads || done)
        {

			Timer totalthread("totalthread", true);
            int numLoops = 0;
            do
            {
				std::cout<<"num_buffers_full : "<<num_buffers_full<<std::endl;
				std::cout<<"done : "<<done<<std::endl;
				
                // Wait for all threads to be ready to receive
                for(int i = 0; i < numThreads; ++i)
                {
					Timer eachthread("eachthread", true);
                    sem_wait(semVec[i]);
                    Thread* pThread = threadVec[i];
                    pThread->swapBuffers(*inputBuffers[i], *outputBuffers[i]);
					std::cout<<"3---------"<<std::endl;
					printf("thread No. %d time : %lf\n", i, eachthread.getElapsedWallTime());
	            }
				
				
                num_buffers_full = 0;
                next_thread = 0;

                // Process the results and clear the buffers
                for(int i = 0; i < numThreads; ++i)
                {
					
					std::cout<<"4---------"<<std::endl;
                    assert(inputBuffers[i]->size() == outputBuffers[i]->size());
                    for(size_t j = 0; j < inputBuffers[i]->size(); ++j)
                    {
						if(i==0 && j==0)
							printf("totalthread time : %lf\n", totalthread.getElapsedWallTime());

						std::cout<<"5---------"<<std::endl;
                        pPostProcessor->process((*inputBuffers[i])[j], (*outputBuffers[i])[j]);
                        ++numWorkItemsWrote;
						std::cout<<"6---------"<<std::endl;
                    }
                    
                    inputBuffers[i]->clear();
                    outputBuffers[i]->clear();
                }
			
                double proc_time_secs = timer.getElapsedWallTime();
                if(generator.getNumConsumed() % (10 * BUFFER_SIZE * numThreads) == 0)
                    printf("[sga] Processed %zu sequences in %lfs (%lf sequences/s)\n", generator.getNumConsumed(), proc_time_secs, (double)generator.getNumConsumed() / proc_time_secs);

                // This should never loop more than twice
                assert(numLoops < 2);
                ++numLoops;
            } while(done && numWorkItemsWrote < numWorkItemsRead);
        }
    }

    // Cleanup
    for(int i = 0; i < numThreads; ++i)
    {
        threadVec[i]->stop(); // Blocks until the thread joins
        delete threadVec[i];

        sem_destroy(semVec[i]);
        delete semVec[i];

        assert(inputBuffers[i]->empty());
        delete inputBuffers[i];

        assert(outputBuffers[i]->empty());
        delete outputBuffers[i];
    }
    assert(n == (size_t)-1 || generator.getNumConsumed() == n);
    assert(numWorkItemsRead == numWorkItemsWrote);

    double proc_time_secs = timer.getElapsedWallTime();
    printf("[sga::process] processed %zu sequences in %lfs (%lf sequences/s)\n", 
            generator.getNumConsumed(), proc_time_secs, (double)generator.getNumConsumed() / proc_time_secs);
    return generator.getNumConsumed();
}

// Design:
// This function is a generic function to read some INPUT from a 
// generic generator object, then perform work on them.
// The actual processing is done by the Processor class 
// that is passed in. The number of threads
// created is determined by the size of the vector of processors - 
// one thread per processor. 
//
// The function buffers batches of input data.
// Once the buffers are full, the reads are dispatched to the thread
// which run the actual processing independently. An optional post processor
// can be specified to process the results that the threads return. If the n
// parameter is used, at most n sequences will be read from the file.
//
// This version is based on OpenMP.
template<class Input, class Output, class Generator, class Processor, class PostProcessor>
size_t processWorkParallelOpenMP(Generator& generator, 
                                 std::vector<Processor*> processPtrVector, 
                                 PostProcessor* pPostProcessor, 
                                 size_t n = -1)
{
#if HAVE_OPENMP
    Timer timer("SequenceProcess", true);

    // Helpful typedefs
    typedef std::vector<Input> InputVector;
    typedef std::vector<Output> OutputVector;

    InputVector inputBuffer;
    OutputVector outputBuffer;

    size_t numWorkItemsRead = 0;
    size_t numWorkItemsWrote = 0;
    size_t numThreads = processPtrVector.size();

    omp_set_num_threads(numThreads);

    bool done = false;
    while(!done)
    {
        // Parse reads from the stream and add them into the incoming buffers
        Input workItem;
        bool valid = generator.generate(workItem);
        if(valid)
        {
            inputBuffer.push_back(workItem);
            numWorkItemsRead += 1;
        }
        
        done = !valid || generator.getNumConsumed() == n;

        // Once all buffers are full or the input is finished, dispatch the work to the threads
        if(inputBuffer.size() == (10 * numThreads * BUFFER_SIZE) || done)
        {
            outputBuffer.resize(inputBuffer.size());

            //
            #pragma omp parallel for schedule(dynamic, 256)
            for(int i = 0; i < (int)inputBuffer.size(); ++i)
            {
                // Dispatch the work to a processor and write the output to the output buffer
                size_t tid = omp_get_thread_num();
                outputBuffer[i] = processPtrVector[tid]->process(inputBuffer[i]);
            }

            // Process the output with a single thread
            for(size_t i = 0; i < inputBuffer.size(); ++i)
            {
                pPostProcessor->process(inputBuffer[i], outputBuffer[i]);
                numWorkItemsWrote += 1;
            }
            inputBuffer.clear();
            outputBuffer.clear();

            double proc_time_secs = timer.getElapsedWallTime();
            if(generator.getNumConsumed() % (numThreads * BUFFER_SIZE) == 0)
                printf("[sga] Processed %zu sequences in %lfs (%lf sequences/s)\n", generator.getNumConsumed(), proc_time_secs, (double)generator.getNumConsumed() / proc_time_secs);
        }
    }

    assert(n == (size_t)-1 || generator.getNumConsumed() == n);
    assert(numWorkItemsRead == numWorkItemsWrote);

    double proc_time_secs = timer.getElapsedWallTime();
    printf("[sga::process] processed %zu sequences in %lfs (%lf sequences/s)\n", 
            generator.getNumConsumed(), proc_time_secs, (double)generator.getNumConsumed() / proc_time_secs);
    return generator.getNumConsumed();
#else // OPENMP
    (void)generator;
    (void)processPtrVector;
    (void)pPostProcessor;
    (void)n;
    printf("Error: threading enabled but you did not compile with OpenMP\n");
    exit(EXIT_FAILURE);
#endif
}

// Wrapper function for operating over n elements of from a SeqReader
template<class Input, class Output, class Processor, class PostProcessor>
size_t processSequencesParallel(SeqReader& reader, 
                                std::vector<Processor*> processPtrVector, 
                                PostProcessor* pPostProcessor, 
                                size_t n = -1)
{
    typedef WorkItemGenerator<Input> InputGenerator;
    InputGenerator generator(&reader);
    return processWorkParallelPthread<Input, 
                                      Output, 
                                      InputGenerator, 
                                      Processor, 
                                      PostProcessor>(generator, processPtrVector, pPostProcessor, n);
}

// Wrapper function for operating over n elements of from a SeqReader
template<class Input, class Output, class Processor, class PostProcessor>
size_t processSequencesParallelOpenMP(SeqReader& reader, 
                                      std::vector<Processor*> processPtrVector, 
                                      PostProcessor* pPostProcessor, 
                                      size_t n = -1)
{
    typedef WorkItemGenerator<Input> InputGenerator;
    InputGenerator generator(&reader);
    return processWorkParallelOpenMP<Input, 
                                     Output, 
                                     InputGenerator, 
                                     Processor, 
                                     PostProcessor>(generator, processPtrVector, pPostProcessor, n);
}

template<class Input, class Output, class Generator, class Processor, class PostProcessor>
size_t processWorkParallelOpenMPbyYi_De(Generator& generator, 
                                 Processor* processPtr, 
                                 PostProcessor* postProcessor, 
                                 size_t n = -1)
{
	
	/*Timer timer("SequenceProcess", true);

    // Helpful typedefs
    typedef std::vector<Input> InputVector;
    typedef std::vector<Output> OutputVector;
	
	Input input;
	Output output;
    InputVector inputBuffer;
    OutputVector outputBuffer;

    size_t numWorkItemsRead = 0;
    size_t numWorkItemsWrote = 0;
    size_t numThreads = processPtrVector.size();

    omp_set_num_threads(numThreads);

    bool done = false;
	size_t i = 0;
    while(i < 1)
    {
        // Parse reads from the stream and add them into the incoming buffers
        Input workItem;
        bool valid = generator[0].generate(workItem);
        if(valid)
        {
            //inputBuffer.push_back(workItem);
            numWorkItemsRead += 1;
        }
        
        done = !valid || generator[0].getNumConsumed() == n;
        output = processPtrVector[0]->process(workItem);
		
		std::cout<<"valid : "<<valid<<std::endl;
		std::cout<<"numWorkItemsRead : "<<numWorkItemsRead<<std::endl;
		//std::cout<<"read : "<<workItem.read.seq.toString()<<std::endl;
		if(done == true)
			std::cout<<"done!!!"<<std::endl;
	
	
        // Process the output with a single thread
        //postProcessorVecor[0]->process(workItem, output);
        numWorkItemsWrote += 1;
   
        /*double proc_time_secs = timer.getElapsedWallTime();
        if(generator[0].getNumConsumed() % (numThreads * BUFFER_SIZE) == 0)
            printf("[sga] Processed %zu sequences in %lfs (%lf sequences/s)\n", generator[0].getNumConsumed(), proc_time_secs, (double)generator[0].getNumConsumed() / proc_time_secs);*/
	//	i++;
        
    //}

    //assert(n == (size_t)-1 || generator[0].getNumConsumed() == n);
    //assert(numWorkItemsRead == numWorkItemsWrote);

    /*double proc_time_secs = timer.getElapsedWallTime();
    printf("[sga::process] processed %zu sequences in %lfs (%lf sequences/s)\n", 
            generator[0].getNumConsumed(), proc_time_secs, (double)generator[0].getNumConsumed() / proc_time_secs);*/
	/*std::cout<<"generator[0].getNumConsumed() : "<<generator[0].getNumConsumed()<<std::endl;*/
    //return generator[0].getNumConsumed();
	
	
	
	Input input;
	Output output;
	size_t numWorkItemsRead = 0;
	bool done =false;
	while(!done)
	{	
	    bool valid = generator.generate(input);
		if(valid)
			numWorkItemsRead += 1;
		done = !valid;
		output = processPtr->process(input);
		postProcessor->process(input, output);
	}
	return generator.getNumConsumed();

}

template<class Input, class Output, class Processor, class PostProcessor>
size_t processSequencesParallelOpenMPbyYi_De(SeqReader &reader, 
                                      Processor* processPtr,
                                      PostProcessor* postProcessor, 
                                      size_t n = -1)
{
	
    typedef WorkItemGenerator<Input> InputGenerator;
	InputGenerator generator(&reader);
	
		
    return processWorkParallelOpenMPbyYi_De<Input, 
                                     Output, 
                                     InputGenerator, 
                                     Processor, 
                                     PostProcessor>(generator, processPtr, postProcessor, n);
	
}


// Wrapper function for operating over a file of sequences
template<class Input, class Output, class Processor, class PostProcessor>
size_t processSequencesParallel(const std::string& readsFile, std::vector<Processor*> processPtrVector, PostProcessor* pPostProcessor)
{
    SeqReader reader(readsFile);
    return processSequencesParallel<Input, Output, Processor, PostProcessor>(reader, processPtrVector, pPostProcessor);
}

// Wrapper function for operating over a file of sequences
template<class Input, class Output, class Processor, class PostProcessor>
size_t processSequencesParallelOpenMP(const std::string& readsFile, std::vector<Processor*> processPtrVector, PostProcessor* pPostProcessor)
{
    SeqReader reader(readsFile);
    return processSequencesParallelOpenMP<Input, Output, Processor, PostProcessor>(reader, processPtrVector, pPostProcessor);
}


template<class Input, class Output, class Processor, class PostProcessor>
size_t processSequencesParallelOpenMPbyYi_De(const std::string& readsFile, Processor* processPtr, PostProcessor* pPostProcessor)
{

	SeqReader reader(readsFile);

    return processSequencesParallelOpenMPbyYi_De<Input, Output, Processor, PostProcessor>(reader, processPtr, pPostProcessor);
}

};

#endif
