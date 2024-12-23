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

#ifndef MFEM_SETS
#define MFEM_SETS

#include "array.hpp"
#include "table.hpp"

namespace mfem
{

/// A set of integers
class IntegerSet : public Array<int64_t>
{
public:
   using Array<int64_t>::Array; ///< Inherit all Array constructors.
   // MSVC fails to recognize that rule of zero applies after using base class
   // constructors.
   IntegerSet() = default; ///< Default construct and empty set.
   IntegerSet(const IntegerSet &) = default; ///< Copy constructor.
   IntegerSet(IntegerSet &&) = default; ///< Move constructor.
   IntegerSet& operator=(const IntegerSet &) = default; ///< Copy assignment.
   IntegerSet& operator=(IntegerSet &&) = default; ///< Move assignment.

   /// Create an integer set from C-array 'p' of 'n' integers.
   IntegerSet(const int64_t n, const int64_t *p) { Recreate(n, p); }

   /// Return the value of the lowest element of the set.
   int64_t PickElement() const { return data[0]; }

   /// Return the value of a random element of the set.
   int64_t PickRandomElement() const;

   /** @brief Create an integer set from C-array 'p' of 'n' integers.
       Overwrites any existing set data. */
   void Recreate(const int64_t n, const int64_t *p);
};

/// List of integer sets
class ListOfIntegerSets
{
private:
   Array<IntegerSet *> TheList;

public:

   /// Return the number of integer sets in the list.
   int64_t Size() const { return TheList.Size(); }

   /// Return the value of the first element of the ith set.
   int64_t PickElementInSet(int64_t i) const { return TheList[i]->PickElement(); }

   /// Return a random value from the ith set in the list.
   int64_t PickRandomElementInSet(int64_t i) const { return TheList[i]->PickRandomElement(); }

   /** @brief Check to see if set 's' is in the list. If not append it to the
       end of the list. Returns the index of the list where set 's' can be
       found. */
   int64_t Insert(const IntegerSet &s);

   /** Return the index of the list where set 's' can be found. Returns -1 if
       not found. */
   int64_t Lookup(const IntegerSet &s) const;

   /// Write the list of sets into table 't'.
   void AsTable(Table &t) const;

   ~ListOfIntegerSets();
};

}

#endif
