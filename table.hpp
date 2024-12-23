// Copyright (c) 2010-2024, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_TABLE
#define MFEM_TABLE

// Data types for Table.

#include "array.hpp"
#include "globals.hpp"
#include <ostream>
#include <istream>

namespace mfem
{

/// Helper struct for defining a connectivity table, see Table::MakeFromList.
struct Connection
{
   int64_t from, to;
   Connection() = default;
   Connection(int64_t from, int64_t to) : from(from), to(to) {}

   bool operator== (const Connection &rhs) const
   { return (from == rhs.from) && (to == rhs.to); }
   bool operator< (const Connection &rhs) const
   { return (from == rhs.from) ? (to < rhs.to) : (from < rhs.from); }
};


/** Data type Table. Table stores the connectivity of elements of TYPE I
    to elements of TYPE II, for example, it may be Element-To-Face
    connectivity table, etc. */
class Table
{
protected:
   /// size is the number of TYPE I elements.
   int64_t size;

   /** Arrays for the connectivity information in the CSR storage.
       I is of size "size+1", J is of size the number of connections
       between TYPE I to TYPE II elements (actually stored I[size]). */
   Memory<int64_t> I, J;

public:
   /// Creates an empty table
   Table() { size = -1; }

   /// Copy constructor
   Table(const Table &);

   /** Merge constructors
       This is used to combine two or three tables into one table.*/
   Table(const Table &table1,
         const Table &table2, int64_t offset2);
   Table(const Table &table1,
         const Table &table2, int64_t offset2,
         const Table &table3, int64_t offset3);

   /// Assignment operator: deep copy
   Table& operator=(const Table &rhs);

   /// Create a table with an upper limit for the number of connections.
   explicit Table (int64_t dim, int64_t connections_per_row = 3);

   /** Create a table from a list of connections, see MakeFromList(). */
   Table(int64_t nrows, Array<Connection> &list) : size(-1)
   { MakeFromList(nrows, list); }

   /** Create a table with one entry per row with column indices given
       by 'partitioning'. */
   Table (int64_t nrows, int64_t *partitioning);

   /// Next 7 methods are used together with the default constructor
   void MakeI (int64_t nrows);
   void AddAColumnInRow (int64_t r) { I[r]++; }
   void AddColumnsInRow (int64_t r, int64_t ncol) { I[r] += ncol; }
   void MakeJ();
   void AddConnection (int64_t r, int64_t c) { J[I[r]++] = c; }
   void AddConnections (int64_t r, const int64_t *c, int64_t nc);
   void ShiftUpI();

   /// Set the size and the number of connections for the table.
   void SetSize(int64_t dim, int64_t connections_per_row);

   /** Set the rows and the number of all connections for the table.
       Does NOT initialize the whole array I ! (I[0]=0 and I[rows]=nnz only) */
   void SetDims(int64_t rows, int64_t nnz);

   /// Returns the number of TYPE I elements.
   inline int64_t Size() const { return size; }

   /** Returns the number of connections in the table. If Finalize() is
       not called, it returns the number of possible connections established
       by the used constructor. Otherwise, it is exactly the number of
       established connections before calling Finalize(). */
   inline int64_t Size_of_connections() const { HostReadI(); return I[size]; }

   /** Returns index of the connection between element i of TYPE I and
       element j of TYPE II. If there is no connection between element i
       and element j established in the table, then the return value is -1. */
   int64_t operator() (int64_t i, int64_t j) const;

   /// Return row i in array row (the Table must be finalized)
   void GetRow(int64_t i, Array<int64_t> &row) const;

   int64_t RowSize(int64_t i) const { return I[i+1]-I[i]; }

   const int64_t *GetRow(int64_t i) const { return J+I[i]; }
   int64_t *GetRow(int64_t i) { return J+I[i]; }

   int64_t *GetI() { return I; }
   int64_t *GetJ() { return J; }
   const int64_t *GetI() const { return I; }
   const int64_t *GetJ() const { return J; }

   Memory<int64_t> &GetIMemory() { return I; }
   Memory<int64_t> &GetJMemory() { return J; }
   const Memory<int64_t> &GetIMemory() const { return I; }
   const Memory<int64_t> &GetJMemory() const { return J; }

   const int64_t *ReadI(bool on_dev = true) const
   { return mfem::Read(I, I.Capacity(), on_dev); }
   int64_t *WriteI(bool on_dev = true)
   { return mfem::Write(I, I.Capacity(), on_dev); }
   int64_t *ReadWriteI(bool on_dev = true)
   { return mfem::ReadWrite(I, I.Capacity(), on_dev); }
   const int64_t *HostReadI() const
   { return mfem::Read(I, I.Capacity(), false); }
   int64_t *HostWriteI()
   { return mfem::Write(I, I.Capacity(), false); }
   int64_t *HostReadWriteI()
   { return mfem::ReadWrite(I, I.Capacity(), false); }

   const int64_t *ReadJ(bool on_dev = true) const
   { return mfem::Read(J, J.Capacity(), on_dev); }
   int64_t *WriteJ(bool on_dev = true)
   { return mfem::Write(J, J.Capacity(), on_dev); }
   int64_t *ReadWriteJ(bool on_dev = true)
   { return mfem::ReadWrite(J, J.Capacity(), on_dev); }
   const int64_t *HostReadJ() const
   { return mfem::Read(J, J.Capacity(), false); }
   int64_t *HostWriteJ()
   { return mfem::Write(J, J.Capacity(), false); }
   int64_t *HostReadWriteJ()
   { return mfem::ReadWrite(J, J.Capacity(), false); }

   /// @brief Sort the column (TYPE II) indices in each row.
   void SortRows();

   /// Replace the #I and #J arrays with the given @a newI and @a newJ arrays.
   /** If @a newsize < 0, then the size of the Table is not modified. */
   void SetIJ(int64_t *newI, int64_t *newJ, int64_t newsize = -1);

   /** Establish connection between element i and element j in the table.
       The return value is the index of the connection. It returns -1 if it
       fails to establish the connection. Possibilities are there is not
       enough memory on row i to establish connection to j, an attempt to
       establish new connection after calling Finalize(). */
   int64_t Push( int64_t i, int64_t j );

   /** Finalize the table initialization. The function may be called
       only once, after the table has been initialized, in order to compress
       array J (by getting rid of -1's in array J). Calling this function
       will "freeze" the table and function Push will work no more.
       Note: The table is functional even without calling Finalize(). */
   void Finalize();

   /** Create the table from a list of connections {(from, to)}, where 'from'
       is a TYPE I index and 'to' is a TYPE II index. The list is assumed to be
       sorted and free of duplicities, i.e., you need to call Array::Sort and
       Array::Unique before calling this method. */
   void MakeFromList(int64_t nrows, const Array<Connection> &list);

   /// Returns the number of TYPE II elements (after Finalize() is called).
   int64_t Width() const;

   /// Call this if data has been stolen.
   void LoseData() { size = -1; I.Reset(); J.Reset(); }

   /// Prints the table to stream out.
   void Print(std::ostream & out = mfem::out, int64_t width = 4) const;
   void PrintMatlab(std::ostream & out) const;

   void Save(std::ostream &out) const;
   void Load(std::istream &in);

   void Copy(Table & copy) const;
   void Swap(Table & other);

   void Clear();

   std::size_t MemoryUsage() const;

   /// Destroys Table.
   ~Table();
};

/// Specialization of the template function Swap<> for class Table
template <> inline void Swap<Table>(Table &a, Table &b)
{
   a.Swap(b);
}

///  Transpose a Table
void Transpose (const Table &A, Table &At, int64_t ncols_A_ = -1);
Table * Transpose (const Table &A);

///  @brief Transpose an Array<int64_t>.
///
/// The array @a A represents a table where each row @a i has exactly one
/// connection to the column (TYPE II) index specified by @a A[i].
///
/// @note The column (TYPE II) indices in each row of @a At will be sorted.
void Transpose(const Array<int64_t> &A, Table &At, int64_t ncols_A_ = -1);

///  C = A * B  (as boolean matrices)
void Mult (const Table &A, const Table &B, Table &C);
Table * Mult (const Table &A, const Table &B);


/** Data type STable. STable is similar to Table, but it's for symmetric
    connectivity, i.e. TYPE I is equivalent to TYPE II. In the first
    dimension we put the elements with smaller index. */
class STable : public Table
{
public:
   /// Creates table with fixed number of connections.
   STable (int64_t dim, int64_t connections_per_row = 3);

   /** Returns index of the connection between element i of TYPE I and
       element j of TYPE II. If there is no connection between element i
       and element j established in the table, then the return value is -1. */
   int64_t operator() (int64_t i, int64_t j) const;

   /** Establish connection between element i and element j in the table.
       The return value is the index of the connection. It returns -1 if it
       fails to establish the connection. Possibilities are there is not
       enough memory on row i to establish connection to j, an attempt to
       establish new connection after calling Finalize(). */
   int64_t Push( int64_t i, int64_t j );

   /// Destroys STable.
   ~STable() {}
};


class DSTable
{
private:
   class Node
   {
   public:
      Node *Prev;
      int64_t  Column, Index;
   };

   int64_t  NumRows, NumEntries;
   Node **Rows;
#ifdef MFEM_USE_MEMALLOC
   MemAlloc <Node, 1024> NodesMem;
#endif

   int64_t Push_(int64_t r, int64_t c);
   int64_t Index(int64_t r, int64_t c) const;

public:
   DSTable(int64_t nrows);
   int64_t NumberOfRows() const { return (NumRows); }
   int64_t NumberOfEntries() const { return (NumEntries); }
   int64_t Push(int64_t a, int64_t b)
   { return ((a <= b) ? Push_(a, b) : Push_(b, a)); }
   int64_t operator()(int64_t a, int64_t b) const
   { return ((a <= b) ? Index(a, b) : Index(b, a)); }
   ~DSTable();

   class RowIterator
   {
   private:
      Node *n;
   public:
      RowIterator (const DSTable &t, int64_t r) { n = t.Rows[r]; }
      int64_t operator!() { return (n != NULL); }
      void operator++() { n = n->Prev; }
      int64_t Column() { return (n->Column); }
      int64_t Index() { return (n->Index); }
      void SetIndex(int64_t new_idx) { n->Index = new_idx; }
   };
};

}

#endif
