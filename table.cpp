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

// Implementation of data types Table.

#include "array.hpp"
#include "table.hpp"
#include "error.hpp"

#include "mem_manager.hpp"
#include <iostream>
#include <iomanip>

namespace mfem
{

using namespace std;

Table::Table(const Table &table)
{
   size = table.size;
   if (size >= 0)
   {
      const int64_t nnz = table.I[size];
      I.New(size+1, table.I.GetMemoryType());
      J.New(nnz, table.J.GetMemoryType());
      I.CopyFrom(table.I, size+1);
      J.CopyFrom(table.J, nnz);
   }
}

Table::Table(const Table &table1,
             const Table &table2, int64_t offset)
{
   MFEM_ASSERT(table1.size == table2.size,
               "Tables have different sizes can not merge.");
   size = table1.size;

   const int64_t nnz = table1.I[size] + table2.I[size];
   I.New(size+1, table1.I.GetMemoryType());
   J.New(nnz, table1.J.GetMemoryType());

   I[0] = 0;
   Array<int64_t> row;
   for (int64_t i = 0; i < size; i++)
   {
      I[i+1] = I[i];

      table1.GetRow(i, row);
      for (int64_t r = 0; r < row.Size(); r++,  I[i+1] ++)
      {
         J[ I[i+1] ] = row[r];
      }

      table2.GetRow(i, row);
      for (int64_t r = 0; r < row.Size(); r++,  I[i+1] ++)
      {
         J[ I[i+1] ] = (row[r] < 0) ? row[r] - offset : row[r] + offset;
      }

   }
}

Table::Table(const Table &table1,
             const Table &table2, int64_t offset2,
             const Table &table3, int64_t offset3)
{
   MFEM_ASSERT(table1.size == table2.size,
               "Tables have different sizes can not merge.");
   MFEM_ASSERT(table1.size == table3.size,
               "Tables have different sizes can not merge.");
   size = table1.size;

   const int64_t nnz = table1.I[size] + table2.I[size] + table3.I[size];
   I.New(size+1, table1.I.GetMemoryType());
   J.New(nnz, table1.J.GetMemoryType());

   I[0] = 0;
   Array<int64_t> row;
   for (int64_t i = 0; i < size; i++)
   {
      I[i+1] = I[i];

      table1.GetRow(i, row);
      for (int64_t r = 0; r < row.Size(); r++,  I[i+1] ++)
      {
         J[ I[i+1] ] = row[r];
      }

      table2.GetRow(i, row);
      for (int64_t r = 0; r < row.Size(); r++,  I[i+1] ++)
      {
         J[ I[i+1] ] = (row[r] < 0) ? row[r] - offset2 : row[r] + offset2;
      }

      table3.GetRow(i, row);
      for (int64_t r = 0; r < row.Size(); r++,  I[i+1] ++)
      {
         J[ I[i+1] ] = (row[r] < 0) ? row[r] - offset3 : row[r] + offset3;
      }
   }
}

Table& Table::operator=(const Table &rhs)
{
   Clear();

   Table copy(rhs);
   Swap(copy);

   return *this;
}

Table::Table (int64_t dim, int64_t connections_per_row)
{
   int64_t i, j, sum = dim * connections_per_row;

   size = dim;
   I.New(size+1);
   J.New(sum);

   I[0] = 0;
   for (i = 1; i <= size; i++)
   {
      I[i] = I[i-1] + connections_per_row;
      for (j = I[i-1]; j < I[i]; j++) { J[j] = -1; }
   }
}

Table::Table (int64_t nrows, int64_t *partitioning)
{
   size = nrows;

   I.New(size+1);
   J.New(size);

   for (int64_t i = 0; i < size; i++)
   {
      I[i] = i;
      J[i] = partitioning[i];
   }
   I[size] = size;
}

void Table::MakeI (int64_t nrows)
{
   SetDims (nrows, 0);

   for (int64_t i = 0; i <= nrows; i++)
   {
      I[i] = 0;
   }
}

void Table::MakeJ()
{
   int64_t i, j, k;

   for (k = i = 0; i < size; i++)
   {
      j = I[i], I[i] = k, k += j;
   }

   J.Delete();
   J.New(I[size]=k);
}

void Table::AddConnections (int64_t r, const int64_t *c, int64_t nc)
{
   int64_t *jp = J+I[r];

   for (int64_t i = 0; i < nc; i++)
   {
      jp[i] = c[i];
   }
   I[r] += nc;
}

void Table::ShiftUpI()
{
   for (int64_t i = size; i > 0; i--)
   {
      I[i] = I[i-1];
   }
   I[0] = 0;
}

void Table::SetSize(int64_t dim, int64_t connections_per_row)
{
   SetDims (dim, dim * connections_per_row);

   if (size > 0)
   {
      I[0] = 0;
      for (int64_t i = 0, j = 0; i < size; i++)
      {
         int64_t end = I[i] + connections_per_row;
         I[i+1] = end;
         for ( ; j < end; j++) { J[j] = -1; }
      }
   }
}

void Table::SetDims(int64_t rows, int64_t nnz)
{
   int64_t j;

   j = (I) ? (I[size]) : (0);
   if (size != rows)
   {
      size = rows;
      I.Delete();
      (rows >= 0) ? I.New(rows+1) : I.Reset();
   }

   if (j != nnz)
   {
      J.Delete();
      (nnz > 0) ? J.New(nnz) : J.Reset();
   }

   if (size >= 0)
   {
      I[0] = 0;
      I[size] = nnz;
   }
}

int64_t Table::operator() (int64_t i, int64_t j) const
{
   if ( i>=size || i<0 )
   {
      return -1;
   }

   int64_t k, end = I[i+1];
   for (k = I[i]; k < end; k++)
   {
      if (J[k] == j)
      {
         return k;
      }
      else if (J[k] == -1)
      {
         return -1;
      }
   }
   return -1;
}

void Table::GetRow(int64_t i, Array<int64_t> &row) const
{
   MFEM_ASSERT(i >= 0 && i < size, "Row index " << i << " is out of range [0,"
               << size << ')');

   HostReadJ();
   HostReadI();

   row.SetSize(RowSize(i));
   row.Assign(GetRow(i));
}

void Table::SortRows()
{
   for (int64_t r = 0; r < size; r++)
   {
      std::sort(J + I[r], J + I[r+1]);
   }
}

void Table::SetIJ(int64_t *newI, int64_t *newJ, int64_t newsize)
{
   I.Delete();
   J.Delete();
   if (newsize >= 0)
   {
      size = newsize;
   }
   I.Wrap(newI, size+1, true);
   J.Wrap(newJ, I[size], true);
}

int64_t Table::Push(int64_t i, int64_t j)
{
   MFEM_ASSERT(i >=0 &&
               i<size, "Index out of bounds.  i = " << i << " size " << size);

   for (int64_t k = I[i], end = I[i+1]; k < end; k++)
   {
      if (J[k] == j)
      {
         return k;
      }
      else if (J[k] == -1)
      {
         J[k] = j;
         return k;
      }
   }

   MFEM_ABORT("Reached end of loop unexpectedly: (i,j) = (" << i << ", " << j
              << ")");

   return -1;
}

void Table::Finalize()
{
   int64_t i, j, end, sum = 0, n = 0, newI = 0;

   for (i=0; i<I[size]; i++)
   {
      if (J[i] != -1)
      {
         sum++;
      }
   }

   if (sum != I[size])
   {
      int64_t *NewJ = Memory<int64_t>(sum);

      for (i=0; i<size; i++)
      {
         end = I[i+1];
         for (j=I[i]; j<end; j++)
         {
            if (J[j] == -1) { break; }
            NewJ[ n++ ] = J[j];
         }
         I[i] = newI;
         newI = n;
      }
      I[size] = sum;

      J.Delete();

      J.Wrap(NewJ, sum, true);

      MFEM_ASSERT(sum == n, "sum = " << sum << ", n = " << n);
   }
}

void Table::MakeFromList(int64_t nrows, const Array<Connection> &list)
{
   Clear();

   size = nrows;
   int64_t nnz = list.Size();

   I.New(size+1);
   J.New(nnz);

   for (int64_t i = 0, k = 0; i <= size; i++)
   {
      I[i] = k;
      while (k < nnz && list[k].from == i)
      {
         J[k] = list[k].to;
         k++;
      }
   }
}

int64_t Table::Width() const
{
   int64_t width = -1, nnz = (size >= 0) ? I[size] : 0;
   for (int64_t k = 0; k < nnz; k++)
   {
      if (J[k] > width) { width = J[k]; }
   }
   return width + 1;
}

void Table::Print(std::ostream & os, int64_t width) const
{
   int64_t i, j;

   for (i = 0; i < size; i++)
   {
      os << "[row " << i << "]\n";
      for (j = I[i]; j < I[i+1]; j++)
      {
         os << setw(5) << J[j];
         if ( !((j+1-I[i]) % width) )
         {
            os << '\n';
         }
      }
      if ((j-I[i]) % width)
      {
         os << '\n';
      }
   }
}

void Table::PrintMatlab(std::ostream & os) const
{
   int64_t i, j;

   for (i = 0; i < size; i++)
   {
      for (j = I[i]; j < I[i+1]; j++)
      {
         os << i << " " << J[j] << " 1. \n";
      }
   }

   os << flush;
}

void Table::Save(std::ostream &os) const
{
   os << size << '\n';

   for (int64_t i = 0; i <= size; i++)
   {
      os << I[i] << '\n';
   }
   for (int64_t i = 0, nnz = I[size]; i < nnz; i++)
   {
      os << J[i] << '\n';
   }
}

void Table::Load(std::istream &in)
{
   I.Delete();
   J.Delete();

   in >> size;
   I.New(size+1);
   for (int64_t i = 0; i <= size; i++)
   {
      in >> I[i];
   }
   int64_t nnz = I[size];
   J.New(nnz);
   for (int64_t j = 0; j < nnz; j++)
   {
      in >> J[j];
   }
}

void Table::Clear()
{
   I.Delete();
   J.Delete();
   size = -1;
   I.Reset();
   J.Reset();
}

void Table::Copy(Table & copy) const
{
   copy = *this;
}

void Table::Swap(Table & other)
{
   mfem::Swap(size, other.size);
   mfem::Swap(I, other.I);
   mfem::Swap(J, other.J);
}

std::size_t Table::MemoryUsage() const
{
   if (size < 0 || I == NULL) { return 0; }
   return (size+1 + I[size]) * sizeof(int64_t);
}

Table::~Table ()
{
   I.Delete();
   J.Delete();
}

void Transpose (const Table &A, Table &At, int64_t ncols_A_)
{
   const int64_t *i_A     = A.GetI();
   const int64_t *j_A     = A.GetJ();
   const int64_t  nrows_A = A.Size();
   const int64_t  ncols_A = (ncols_A_ < 0) ? A.Width() : ncols_A_;
   const int64_t  nnz_A   = i_A[nrows_A];

   At.SetDims (ncols_A, nnz_A);

   int64_t *i_At = At.GetI();
   int64_t *j_At = At.GetJ();

   for (int64_t i = 0; i <= ncols_A; i++)
   {
      i_At[i] = 0;
   }
   for (int64_t i = 0; i < nnz_A; i++)
   {
      i_At[j_A[i]+1]++;
   }
   for (int64_t i = 1; i < ncols_A; i++)
   {
      i_At[i+1] += i_At[i];
   }

   for (int64_t i = 0; i < nrows_A; i++)
   {
      for (int64_t j = i_A[i]; j < i_A[i+1]; j++)
      {
         j_At[i_At[j_A[j]]++] = i;
      }
   }
   for (int64_t i = ncols_A; i > 0; i--)
   {
      i_At[i] = i_At[i-1];
   }
   i_At[0] = 0;
}

Table * Transpose(const Table &A)
{
   Table * At = new Table;
   Transpose(A, *At);
   return At;
}

void Transpose(const Array<int64_t> &A, Table &At, int64_t ncols_A_)
{
   At.MakeI((ncols_A_ < 0) ? (A.Max() + 1) : ncols_A_);
   for (int64_t i = 0; i < A.Size(); i++)
   {
      At.AddAColumnInRow(A[i]);
   }
   At.MakeJ();
   for (int64_t i = 0; i < A.Size(); i++)
   {
      At.AddConnection(A[i], i);
   }
   At.ShiftUpI();
}

void Mult (const Table &A, const Table &B, Table &C)
{
   int64_t  i, j, k, l, m;
   const int64_t *i_A     = A.GetI();
   const int64_t *j_A     = A.GetJ();
   const int64_t *i_B     = B.GetI();
   const int64_t *j_B     = B.GetJ();
   const int64_t  nrows_A = A.Size();
   const int64_t  nrows_B = B.Size();
   const int64_t  ncols_A = A.Width();
   const int64_t  ncols_B = B.Width();

   MFEM_VERIFY( ncols_A <= nrows_B, "Table size mismatch: ncols_A = " << ncols_A
                << ", nrows_B = " << nrows_B);

   Array<int64_t> B_marker (ncols_B);

   for (i = 0; i < ncols_B; i++)
   {
      B_marker[i] = -1;
   }

   int64_t counter = 0;
   for (i = 0; i < nrows_A; i++)
   {
      for (j = i_A[i]; j < i_A[i+1]; j++)
      {
         k = j_A[j];
         for (l = i_B[k]; l < i_B[k+1]; l++)
         {
            m = j_B[l];
            if (B_marker[m] != i)
            {
               B_marker[m] = i;
               counter++;
            }
         }
      }
   }

   C.SetDims (nrows_A, counter);

   for (i = 0; i < ncols_B; i++)
   {
      B_marker[i] = -1;
   }

   int64_t *i_C = C.GetI();
   int64_t *j_C = C.GetJ();
   counter = 0;
   for (i = 0; i < nrows_A; i++)
   {
      i_C[i] = counter;
      for (j = i_A[i]; j < i_A[i+1]; j++)
      {
         k = j_A[j];
         for (l = i_B[k]; l < i_B[k+1]; l++)
         {
            m = j_B[l];
            if (B_marker[m] != i)
            {
               B_marker[m] = i;
               j_C[counter] = m;
               counter++;
            }
         }
      }
   }
}


Table * Mult (const Table &A, const Table &B)
{
   Table * C = new Table;
   Mult(A,B,*C);
   return C;
}

STable::STable (int64_t dim, int64_t connections_per_row) :
   Table(dim, connections_per_row)
{}

int64_t STable::operator() (int64_t i, int64_t j) const
{
   if (i < j)
   {
      return Table::operator()(i,j);
   }
   else
   {
      return Table::operator()(j,i);
   }
}

int64_t STable::Push( int64_t i, int64_t j )
{
   if (i < j)
   {
      return Table::Push(i, j);
   }
   else
   {
      return Table::Push(j, i);
   }
}


DSTable::DSTable(int64_t nrows)
{
   Rows = new Node*[nrows];
   for (int64_t i = 0; i < nrows; i++)
   {
      Rows[i] = NULL;
   }
   NumRows = nrows;
   NumEntries = 0;
}

int64_t DSTable::Push_(int64_t r, int64_t c)
{
   MFEM_ASSERT(r >= 0 && r < NumRows,
               "Row out of bounds: r = " << r << ", NumRows = " << NumRows);
   Node *n;
   for (n = Rows[r]; n != NULL; n = n->Prev)
   {
      if (n->Column == c)
      {
         return (n->Index);
      }
   }
#ifdef MFEM_USE_MEMALLOC
   n = NodesMem.Alloc ();
#else
   n = new Node;
#endif
   n->Column = c;
   n->Index  = NumEntries;
   n->Prev   = Rows[r];
   Rows[r]   = n;
   return (NumEntries++);
}

int64_t DSTable::Index(int64_t r, int64_t c) const
{
   MFEM_ASSERT( r>=0, "Row index must be non-negative, not "<<r);
   if (r >= NumRows)
   {
      return (-1);
   }
   for (Node *n = Rows[r]; n != NULL; n = n->Prev)
   {
      if (n->Column == c)
      {
         return (n->Index);
      }
   }
   return (-1);
}

DSTable::~DSTable()
{
#ifdef MFEM_USE_MEMALLOC
   // NodesMem.Clear();  // this is done implicitly
#else
   for (int64_t i = 0; i < NumRows; i++)
   {
      Node *na, *nb = Rows[i];
      while (nb != NULL)
      {
         na = nb;
         nb = nb->Prev;
         delete na;
      }
   }
#endif
   delete [] Rows;
}

}
