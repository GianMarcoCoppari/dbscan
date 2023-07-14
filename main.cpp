#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <array>
#include <algorithm>
#include <random>


// operator function for necessary operation between data
// only declared here, implemented after main function
template <typename T, std::size_t N> std::array<T, N>& operator -= (std::array<T, N>&, std::array<T, N> const&);
template <typename T, std::size_t N> std::array<T, N> operator - (std::array<T, N> const&, std::array<T, N> const&);
template <typename T, std::size_t N> T operator * (std::array<T, N> const&, std::array<T, N> const&);
template <typename T, std::size_t N, typename K> std::array<T, N>& operator *= (std::array<T, N>&, K const&);
template <typename T, std::size_t N, typename K> std::array<T, N> operator * (K const&, std::array<T, N> const&);
template <typename T, std::size_t N> std::ostream& operator << (std::ostream&, std::array<T, N> const&);


template <typename T, std::size_t N>
std::vector<std::array<T, N>> ReadCSV(std::string const& filepath) {
	/// legge un file csv che contiene dati numerici
	/// il numero di features deve essere fornito dall'utente
	/// e deve essere un numero valido a seconda del 
	/// file di testo fornito. 
	/// Nessun controllo viene effettuato, quindi si assume che l'utente non scriva porcherie
	/// 
	std::fstream fin{filepath.c_str(), std::ios_base::in}; // apro il file in lettura
	
	std::string line{};					// string che contiene la riga del file di testo.
	std::string word{};					// variabile che contiene temporaneamente una parola della stringa precedente
	std::vector<std::string> words{};	// vettore che contiene le parole
	
	std::array<T, N> entry{};					// variabile per gestire una riga del database
	std::vector<std::array<T, N>> entries{};	// dataset finale
	
	while (!fin.eof()) {									// finché non arrivo alla fine del file
		words.clear();					// setup delle variabili all'inizio di ogni ciclo
		entry = std::array<T, N> {};
		
		std::getline(fin, line, '\n');	// leggo una riga dal file 'fin' e la salvo nella variabile 'line', finché non trovo il carattere '\n'.
		std::stringstream str(line);	// creo un stream per la lettura dei dati (parole)
		
		while (std::getline(str, word, ',')) {
			// dallo stream precedente leggo i caratteri
			// finché non incontro il carattere ','
			// ogni parola letta viene inserita quindi nel vettoew 'words'
			words.push_back(word);
		}
		
		for (std::size_t fcount{}; fcount < N; fcount++) {
			// ogni parola viene convertita in formato numerico per essere gestita dal database
			entry[fcount] = std::stod(words[fcount]);
		}
		
		// inserisco la variabile nel dataset
		entries.push_back(entry);
	}
	
	return entries;
}


// toy-dataset generator function to make the dbscan algorithm work
// MakeBolobs generates normally distributed blobs at different centers
// MakeCircleCrown generates concentric circular crowns 
//     1. r is normally distributed around the mean value (parameter)
//     2. theta is uniformly distributed
auto MakeBlobs(std::vector<std::array<double, 2U>> const& centers, double const& stddev, std::size_t const& seed, std::size_t const& nsamples) {
	std::vector<std::array<double, 2U>> dataset{};
	
	std::random_device rd{};
	std::seed_seq seq{seed};
	std::mt19937_64 mt{};
	
	for (auto const& center : centers) {
		std::normal_distribution<double> ux{center[0U], stddev};
		std::normal_distribution<double> uy{center[1U], stddev};
	
		for (std::size_t i{}; i != nsamples; i++) {
			dataset.push_back(std::array<double, 2U> {ux(mt), uy(mt)});
		}
		
	}
	
	return dataset;
}
auto MakeCircleCrown(double const& r, double const& R, double const& stddev, std::size_t const& seed, std::size_t const& nsamples) {
	std::random_device rd{};
	std::seed_seq seq{ seed };
	std::mt19937_64 mt{ seq };
	
	std::normal_distribution<double> inner{r, stddev};
	std::normal_distribution<double> outer{R, stddev};
	std::uniform_real_distribution<double> angle{0, 2 * std::acos(-1)};
	double theta;
		
	std::vector<std::array<double, 2U>> dataset{};
	
	for (std::size_t i{}; i != nsamples; i++) {
		theta = angle(mt);
		dataset.push_back(inner(mt) * std::array<double, 2U> {std::cos(theta), std::sin(theta)});
	}
	for (std::size_t i{}; i != nsamples; i++) {
		theta = angle(mt);
		dataset.push_back(outer(mt) * std::array<double, 2U> {std::cos(theta), std::sin(theta)});
	}
	
	return dataset;
}


int main() {
	// Auxiliary functions
	//     1. subset: computes a subset of the vector given as parameter.
	//	              returns a vector containing indices of elements in the original vector 
	//                that fulfill the condition imposed by parameter 'query'.
	//     2. Union: it is a shell for 'std::set_unioin' function from the STL.
	//               return the union set of the containers given as parameters.
	auto subset = [] <typename DataType, typename Query> (DataType const& Q, std::vector<DataType> const& dataset, Query const& query) {
		std::vector<std::size_t> s{};
		for (auto const& entry : dataset) {
			if (query(Q, entry)) {
				s.push_back(&entry - &dataset[0U]);
			}
		}
		
		return s;
	};
	auto Union = [] <typename DataType> (std::vector<DataType> const& A, std::vector<DataType> const& B) {
		std::vector<DataType> U{};
		std::set_union(A.begin(), A.end(), B.begin(), B.end(), std::back_inserter(U));
		
		return U;
	};
	

	auto DBSCAN = [subset, Union] <typename DataType, typename Query> (std::vector<DataType> const& dataset, std::size_t const& minsize, Query const& query) {
		/// DBSCAN Algorithm imlemented here.
		/// template parameter allow the algorithm to adapt to different types of dataset (in number of features)
		/// numerical values are the only accepted ones. 
		/// this guarantees stability with respect to the dimensionality of the feature space
		
		/// template parameter 'Query' can be a function or a lambda-expression to measure distances in feature space
		/// template allows for passing metric functions (or lambda-expression) with different signature without changing the function

		// cluster counter: 
		//     1. conuts how many clusters are detected
		//     2. numerically, give the label for the cluster
		int ccount = 0; 

		int const noise = -1;		// label for 'noise' points in the dataset
		int const undefined = 0;	// default label for unclassified data
		
		std::vector<int> labels{};		// vector containing the label associated to the corresponding point in 'datatset'
		labels.resize(dataset.size());	// default initialization of all labels to 'undefined'
		
		typedef std::vector<DataType> Cluster;	// typedef for code clarity
		std::vector<Cluster> clusters{};		// clusters vector: points belonging to the same cluster
												// are stored in the same object
		
		for (auto const& entry : dataset) {
			auto const& j = &entry - &dataset[0U]; // index of current element
			
			
			if (labels[j] != undefined) { continue; } // already labelled point, skip
			
			auto group = subset(entry, dataset, query);
			if (group.size() < minsize) { // not enough neighbours
				labels[j] = noise;		  // labels this point as 'noise'
				continue;				  // go to next point
			}
			
			ccount++;			// sufficient number of neighbours -> core point. Update cluster counter value
			labels[j] = ccount; // update label of point j
			
			clusters.push_back(Cluster{});			// create an empty Cluster object
			clusters.back().push_back(dataset[j]);	// add the point to the last dataset
			
			// while the neighbour vector is not empty
			while(!group.empty()) {
				auto const& last = group.back();	// index (in 'datatset') of last neighbour
				group.pop_back();					// delete last element, it will be classified anyway
													// there is no more need to keep track of it
				
				// check element's label value
				if (labels[last] == noise) {
					// if noise, update to the current label value and move on
					labels[last] = ccount;
					clusters.back().push_back(dataset[last]);
					continue;
				}
				
				if (labels[last] != undefined) { continue; }	// skip if already classified
				
				// if last neighbour has not been classified yet -> update its label to current value
				labels[last] = ccount;
				clusters.back().push_back(dataset[last]);
				
				// compute neighbours to neighbour point
				auto const& nearneighbours = subset(dataset[last], dataset, query);
				if (nearneighbours.size() >= minsize) {
					// if suffient number of neighbours-to-neighbour -> neighbour is the union set of the two groups 
					group = Union(group, nearneighbours);
				}
			}
		}
		
		typedef std::pair<std::vector<Cluster>, std::vector<int>> ClusteringResult;
		return ClusteringResult{clusters, labels};
	};
	
	// ------------------------------------------------------------------------------------------------------------------------
	// Next two blocks generate toy datasets
	// Decomment the one you prefer, comment the other
	
	// Block 1 - Blob Dataset
	// std::vector<std::array<double, 2U>> const centers {std::array<double, 2U> {1,   1}, std::array<double, 2U> {1,  -1}, std::array<double, 2U> {-1, -1}};
	// double const stddev{ 0.25 };
	// std::size_t const seed{ 0 };
	// std::size_t const nsamples{ 500 };
	// auto const dataset = MakeBlobs(centers, stddev, seed, nsamples);

	// Block 2 - Circle Crown Dataset
	double const r{ 1 };
	double const R{ 5 };
	double const stddev{ 0.25 };
	std::size_t const seed{ 0 };
	std::size_t const nsamples{ 500 };
	auto const dataset = MakeCircleCrown(r, R, stddev, seed, nsamples);
	

	// ---------------------------------------------------------------
	// DBSCAN Algoritm in action.
	// 1. Set hyper-parameters
	double const epsilon{ 0.5 };	// raggio del 'vicinato'
	std::size_t const minsize{ 8 }; // nuemero minimo di elementi per costituire un cluster

	// euclidean metric for computing distances between points
	auto EuclideanDistance2 = [] <typename Point> (Point const& lhs, Point const& rhs) {
		return (lhs - rhs) * (lhs - rhs);
	};
	// query function for DBSCAN Algorithm
	auto IsInEpsilonVicinity = [epsilon, EuclideanDistance2] <typename Point> (Point const& lhs, Point const& rhs) {
		return EuclideanDistance2(lhs, rhs) < epsilon * epsilon;
	};


	// first output file - dataset
	std::fstream fout{"./concentric-circles.data", std::ios_base::out};
	fout << "x\ty" << std::endl;
	for (auto const& entry : dataset) {
		fout << entry[0U] << '\t' << entry[1U] << std::endl;
	}
	
	auto const result = DBSCAN(dataset, minsize, IsInEpsilonVicinity); // invoco la funzione per generare i cluster e le etichette
	auto const& labels = result.second;
	
	// other output files - one file for each detected cluster
	std::fstream file;
	std::string const filename{"./cluster-"};
	auto const& clusters = result.first;
	for (auto const& cluster : clusters) {
		file.open(filename + std::to_string(&cluster - &clusters[0U]) + std::string{".data"}, std::ios_base::out);
		if (!file.is_open()) { continue; }
		
		file << "x\ty" << std::endl;
		for (auto const& element : cluster) {
			file << element[0U] << '\t' << element[1U] << std::endl;
		}
		
		file.clear();
		file.close();
	}


	return 0;
}

// array operator function implementation -------------------------------------------
template <typename T, std::size_t N> std::array<T, N>& operator -= (std::array<T, N>& lhs, std::array<T, N> const& rhs) {
	for (auto& i : lhs) {
		i -= rhs[&i - &lhs[0U]];
	}
	
	return lhs;
}
template <typename T, std::size_t N> std::array<T, N> operator - (std::array<T, N> const& lhs, std::array<T, N> const& rhs) {
	auto res{lhs};
	return res -= rhs;
}
template <typename T, std::size_t N> T operator * (std::array<T, N> const& lhs, std::array<T, N> const& rhs) {
	T sum{};
	for (auto const& i : lhs) {
		sum += i * rhs[&i - &lhs[0U]];
	}
	
	return sum;
}
template <typename T, std::size_t N, typename K> std::array<T, N>& operator *= (std::array<T, N>& arr, K const& k) {
	for (auto& i : arr) {
		i *= k;
	}
	
	return arr;
}
template <typename T, std::size_t N, typename K> std::array<T, N> operator * (K const& k, std::array<T, N> const& arr) {
	auto res{arr};
	return res *= k;
}
template <typename T, std::size_t N> std::ostream& operator << (std::ostream& os, std::array<T, N> const& arr) {
	for (auto const& i : arr) {
		os << i << "  ";
	}
	
	return os;
}
