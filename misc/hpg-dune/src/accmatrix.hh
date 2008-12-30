#ifndef ACCMATRIX_HH
#define ACCMATRIX_HH

#include <vector>
#include <algorithm>

template<typename ct>
class AccumulationMatrix
{
public:
  enum State {
    PREPARE,
    FINISHED
  };

public:
  struct Triple {
    int rowIndex;
    int columnIndex;
    ct value;

    Triple(int rowInd, int colInd, ct val=0.0) : rowIndex(rowInd), columnIndex(colInd), value(val){}

    bool operator< (const Triple& other) const {
      return rowIndex<other.rowIndex || 
	rowIndex==other.rowIndex && columnIndex<other.columnIndex;
    }
  };

  struct ValueAccessor {
    AccumulationMatrix& matrix;
    const int rowIndex, columnIndex;

    ValueAccessor(AccumulationMatrix &m, int rowInd, int colInd) : matrix(m), rowIndex(rowInd), columnIndex(colInd) 
    { }

    void operator=(ct value) {
      if (matrix.state == PREPARE) {
	matrix.add(rowIndex, columnIndex, value);
      }
    }

    void operator+=(ct value) {
      if (matrix.state == PREPARE) {
	matrix.add(rowIndex, columnIndex, value);
      }
    }

    operator double () {
      return matrix.get(rowIndex, columnIndex);
    }
  };

  struct RowAccessor {
    AccumulationMatrix& matrix;
    const int rowIndex;
    RowAccessor(AccumulationMatrix& m, int rowInd) : matrix(m), rowIndex(rowInd) { }

    const ct& operator[](int colInd) const {
      return matrix._vals[colInd];
    }

    ValueAccessor operator[](int colInd) {
      return ValueAccessor(matrix, rowIndex, colInd);
    }
  };

public:
  typedef ct ctype;

  AccumulationMatrix(int nrows, int ncols) : state(PREPARE), _nrows(nrows), _ncols(ncols)
  {
  }

  RowAccessor operator[](int rowInd) {
    return RowAccessor(*this, rowInd);
  }

  void add(int rowInd, int colInd, ct value) {
    _triples.push_back(Triple(rowInd, colInd, value));
  }

  ct get(int rowInd, int colInd) {
    if (state == PREPARE) {
      for (int i=0; i<_triples.size(); i++) {
	Triple &triple = _triples[i];
	if (triple.rowIndex==rowInd && triple.columnIndex==colInd)
	  return triple.value;
      } 
    }else {
      typedef typename std::vector<Triple>::iterator Iterator;
      Triple value(rowInd, colInd);
      std::pair<Iterator,Iterator> answer =
	std::equal_range(_triples.begin(), _triples.end(), value);
      if (answer.first != answer.second)
	return answer.first->value;
    }
    return 0.0;
  }

  typedef typename std::vector<Triple>::iterator TripleIterator;

  /*! sort and shrink as necessary */
  void finish() {
    std::sort(_triples.begin(), _triples.end());
    std::vector<Triple> shrinked;
    typename std::vector<Triple>::iterator endTriples = _triples.end();
    for (typename std::vector<Triple>::iterator it = _triples.begin(); it!= endTriples; ) {
      Triple triple = *it;
      it++;
      while (it!=endTriples && !(triple<*it)) {
	triple.value+=it->value;
	it++;
      } 
      if (std::abs(triple.value) > 10E-10)
	shrinked.push_back(triple);
    }
    _triples = shrinked;
    state = FINISHED;
  }

  TripleIterator beginTriples() { return _triples.begin(); }
  TripleIterator endTriples() { return _triples.end(); }
  int nnz() { return _triples.size(); }

  int rowdim() { return _nrows; }
  int coldim() { return _ncols; }

private:
  State state;
  
  int _nrows, _ncols;
  std::vector<Triple> _triples;

};

#endif
