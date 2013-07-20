#ifndef edwenger_HiVertexAnalyzer_HiVertexComparator_
#define edwenger_HiVertexAnalyzer_HiVertexComparator_h

template<typename T>
struct GreaterByTracksSize {
	typedef T first_argument_type;
	typedef T second_argument_type;
	bool operator()( const T & t1, const T & t2 ) const {
		return t1.tracksSize() > t2.tracksSize();
	}
};

template<typename T>
struct LessByNormalizedChi2 {
	typedef T first_argument_type;
	typedef T second_argument_type;
	bool operator()( const T & t1, const T & t2 ) const {
		return t1.normalizedChi2() < t2.normalizedChi2();
	}
};

template<typename T>
struct MoreTracksThenLowerChi2 {
	typedef T first_argument_type;
	typedef T second_argument_type;
	bool operator()( const T & t1, const T & t2 ) const {
	  if(t1.tracksSize() == t2.tracksSize())
	    return t1.normalizedChi2() < t2.normalizedChi2();
	  else
	    return t1.tracksSize() > t2.tracksSize();
	}
};

#endif
