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

#ifndef MFEM_MESH
#define MFEM_MESH

#include "stable3d.hpp"
#include "globals.hpp"
#include "attribute_sets.hpp"
#include "element.hpp"
#include "vertex.hpp"
#include "vtk.hpp"

#include <iostream>
#include <array>
#include <map>
#include <vector>
#include <memory>

namespace mfem
{

/** An enum type to specify if interior or boundary faces are desired. */
enum class FaceType : bool {Interior, Boundary};

/// Mesh data type
class Mesh
{
protected:
   int Dim;
   int spaceDim;

   int NumOfVertices, NumOfElements, NumOfBdrElements;
   int NumOfEdges, NumOfFaces;
   /** These variables store the number of Interior and Boundary faces. Calling
       fes->GetMesh()->GetNBE() doesn't return the expected value in 3D because
       periodic meshes in 3D have some of their faces marked as boundary for
       visualization purpose in GLVis. */
   mutable int nbInteriorFaces, nbBoundaryFaces;

   // see MeshGenerator(); global in parallel
   int meshgen;
   // sum of (1 << geom) for all geom of all dimensions; local in parallel
   int mesh_geoms;

   Array<Element *> elements;
   // Vertices are only at the corners of elements, where you would expect them
   // in the lowest-order mesh. In some cases, e.g. in a Mesh that defines the
   // patch topology for a NURBS mesh (see LoadPatchTopo()) the vertices may be
   // empty while NumOfVertices is positive.
   Array<Vertex> vertices;
   Array<Element *> boundary;
   Array<Element *> faces;

   /** @brief This structure stores the low level information necessary to
       interpret the configuration of elements on a specific face. This
       information can be accessed using methods like GetFaceElements(),
       GetFaceInfos(), FaceIsInterior(), etc.

       For accessing higher level deciphered information look at
       Mesh::FaceInformation, and its accessor Mesh::GetFaceInformation().

       Each face contains information on the indices, local reference faces,
       orientations, and potential nonconformity for the two neighboring
       elements on a face.
       Each face can either be an interior, boundary, or shared interior face.
       Each interior face is shared by two elements referred as Elem1 and Elem2.
       For boundary faces only the information on Elem1 is relevant.
       Shared interior faces correspond to faces where Elem1 and Elem2 are
       distributed on different MPI ranks.
       Regarding conformity, three cases are distinguished, conforming faces,
       nonconforming slave faces, and nonconforming master faces. Master and
       slave referring to the coarse and fine elements respectively on a
       nonconforming face.
       Nonconforming slave faces always have the slave element as Elem1 and
       the master element as Elem2. On the other side, nonconforming master
       faces always have the master element as Elem1, and one of the slave
       element as Elem2. Except for ghost nonconforming slave faces, where
       Elem1 is the master side and Elem2 is the slave side.

       The indices of Elem1 and Elem2 can be indirectly extracted from
       FaceInfo::Elem1No and FaceInfo::Elem2No, read the note below for special
       cases on the index of Elem2.

       The local face identifiers are deciphered from FaceInfo::Elem1Inf and
       FaceInfo::Elem2Inf through the formula: LocalFaceIndex = ElemInf/64,
       the semantic of the computed local face identifier can be found in
       fem/geom.cpp. The local face identifier corresponds to an index
       in the Constants<Geometry>::Edges arrays for 2D element geometries, and
       to an index in the Constants<Geometry>::FaceVert arrays for 3D element
       geometries.

       The orientation of each element relative to a face is obtained through
       the formula: Orientation = ElemInf%64, the semantic of the orientation
       can also be found in fem/geom.cpp. The orientation corresponds to
       an index in the Constants<Geometry>::Orient arrays, providing the
       sequence of vertices identifying the orientation of an edge/face. By
       convention the orientation of Elem1 is always set to 0, serving as the
       reference orientation. The orientation of Elem2 relatively to Elem1 is
       therefore determined just by using the orientation of Elem2. An important
       special case is the one of nonconforming faces, the orientation should
       be composed with the PointMatrix, which also contains orientation
       information. A special treatment should be done for 2D, the orientation
       in the PointMatrix is not included, therefore when applying the
       PointMatrix transformation, the PointMatrix should be flipped, except for
       shared nonconforming slave faces where the transformation can be applied
       as is.

       Another special case is the case of shared nonconforming faces. Ghost
       faces use a different design based on so called "ghost" faces.
       Ghost faces, as their name suggest are very well hidden, and they
       usually have a separate interface from "standard" faces.
   */
   struct FaceInfo
   {
      // Inf = 64 * LocalFaceIndex + FaceOrientation
      int Elem1No, Elem2No, Elem1Inf, Elem2Inf;
      int NCFace; /* -1 if this is a regular conforming/boundary face;
                     index into 'nc_faces_info' if >= 0. */
   };
   // NOTE: in NC meshes, master faces have Elem2No == -1. Slave faces on the
   // other hand have Elem2No and Elem2Inf set to the master face's element and
   // its local face number.
   //
   // A local face is one generated from a local element and has index i in
   // faces_info such that i < GetNumFaces(). Also, Elem1No always refers to the
   // element (slave or master, in the nonconforming case) that generated the
   // face.
   // Classification of a local (non-ghost) face based on its FaceInfo:
   // - Elem2No >= 0 --> local interior face; can be either:
   //    - NCFace == -1 --> conforming face, or
   //    - NCFace >= 0 --> nonconforming slave face; Elem2No is the index of
   //      the master volume element; Elem2Inf%64 is 0, see the note in
   //      Mesh::GenerateNCFaceInfo().
   // - Elem2No < 0 --> local "boundary" face; can be one of:
   //    - NCFace == -1 --> conforming face; can be either:
   //       - Elem2Inf < 0 --> true boundary face (no element on side 2)
   //       - Elem2Inf >= 0 --> shared face where element 2 is a face-neighbor
   //         element with index -1-Elem2No. This state is initialized by
   //         ParMesh::ExchangeFaceNbrData().
   //    - NCFace >= 0 --> nonconforming face; can be one of:
   //       - Elem2Inf < 0 --> master nonconforming face, interior or shared;
   //         In this case, Elem2No is -1; see GenerateNCFaceInfo().
   //       - Elem2Inf >= 0 --> shared slave nonconforming face where element 2
   //         is the master face-neighbor element with index -1-Elem2No; see
   //         ParNCMesh::GetFaceNeighbors().
   //
   // A ghost face is a nonconforming face that is generated by a non-local,
   // i.e. ghost, element. A ghost face has index i in faces_info such that
   // i >= GetNumFaces().
   // Classification of a ghost (non-local) face based on its FaceInfo:
   // - Elem1No == -1 --> master ghost face? These ghost faces also have:
   //   Elem2No == -1, Elem1Inf == Elem2Inf == -1, and NCFace == -1.
   // - Elem1No >= 0 --> slave ghost face; Elem1No is the index of the local
   //   master side element, i.e. side 1 IS NOT the side that generated the
   //   face. Elem2No is < 0 and -1-Elem2No is the index of the ghost
   //   face-neighbor element that generated this slave ghost face. In this
   //   case, Elem2Inf >= 0 and NCFace >= 0.
   // Relevant methods: GenerateFaces(), GenerateNCFaceInfo(),
   //                   ParNCMesh::GetFaceNeighbors(),
   //                   ParMesh::ExchangeFaceNbrData()

   struct NCFaceInfo
   {
      bool Slave; // true if this is a slave face, false if master face
      int MasterFace; // if Slave, this is the index of the master face
      // If not Slave, 'MasterFace' is the local face index of this master face
      // as a face in the unique adjacent element.
      const DenseMatrix* PointMatrix; // if Slave, position within master face
      // (NOTE: PointMatrix points to a matrix owned by NCMesh.)

      NCFaceInfo() = default;

      NCFaceInfo(bool slave, int master, const DenseMatrix* pm)
         : Slave(slave), MasterFace(master), PointMatrix(pm) {}
   };

   Array<FaceInfo> faces_info;
   Array<NCFaceInfo> nc_faces_info;

   Table *el_to_edge;
   Table *el_to_face;
   Table *el_to_el;
   Array<int> be_to_face; // faces = vertices (1D), edges (2D), faces (3D)

   Table *bel_to_edge;    // for 3D only

   // Note that the following tables are owned by this class and should not be
   // deleted by the caller. Of these three tables, only face_edge and
   // edge_vertex are returned by access functions.
   mutable Table *face_to_elem;  // Used by FindFaceNeighbors, not returned.
   mutable Table *face_edge;     // Returned by GetFaceEdgeTable().
   mutable Table *edge_vertex;   // Returned by GetEdgeVertexTable().

   static const int vtk_quadratic_tet[10];
   static const int vtk_quadratic_pyramid[13];
   static const int vtk_quadratic_wedge[18];
   static const int vtk_quadratic_hex[27];

public:
   typedef Geometry::Constants<Geometry::SEGMENT>     seg_t;
   typedef Geometry::Constants<Geometry::TRIANGLE>    tri_t;
   typedef Geometry::Constants<Geometry::SQUARE>      quad_t;
   typedef Geometry::Constants<Geometry::TETRAHEDRON> tet_t;
   typedef Geometry::Constants<Geometry::CUBE>        hex_t;
   typedef Geometry::Constants<Geometry::PRISM>       pri_t;
   typedef Geometry::Constants<Geometry::PYRAMID>     pyr_t;

   enum Operation { NONE, REFINE, DEREFINE, REBALANCE };

   /// A list of all unique element attributes used by the Mesh.
   Array<int> attributes;
   /// A list of all unique boundary attributes used by the Mesh.
   Array<int> bdr_attributes;

   /// Named sets of element attributes
   AttributeSets attribute_sets;

   /// Named sets of boundary element attributes
   AttributeSets bdr_attribute_sets;

   // Global parameter that can be used to control the removal of unused
   // vertices performed when reading a mesh in MFEM format. The default value
   // (true) is set in mesh_readers.cpp.
   static bool remove_unused_vertices;

protected:
   void Init();
   void InitTables();
   void SetEmpty();  // Init all data members with empty values
   void DestroyTables();
   void DeleteTables() { DestroyTables(); InitTables(); }
   void DestroyPointers(); // Delete data specifically allocated by class Mesh.
   void Destroy();         // Delete all owned data.
   void ResetLazyData();

   Element *ReadElementWithoutAttr(std::istream &input);
   static void PrintElementWithoutAttr(const Element *el, std::ostream &os);

   Element *ReadElement(std::istream &input);
   static void PrintElement(const Element *el, std::ostream &os);

   // Readers for different mesh formats, used in the Load() method.
   // The implementations of these methods are in mesh_readers.cpp.
   void ReadMFEMMesh(std::istream &input, int version, int &curved);
   void ReadLineMesh(std::istream &input);
   void ReadNetgen2DMesh(std::istream &input, int &curved);
   void ReadNetgen3DMesh(std::istream &input);
   void ReadTrueGridMesh(std::istream &input);
   void CreateVTKMesh(const Vector &points, const Array<int> &cell_data,
                      const Array<int> &cell_offsets,
                      const Array<int> &cell_types,
                      const Array<int> &cell_attributes,
                      int &curved, int &read_gf, bool &finalize_topo);
   void ReadVTKMesh(std::istream &input, int &curved, int &read_gf,
                    bool &finalize_topo);
   void ReadXML_VTKMesh(std::istream &input, int &curved, int &read_gf,
                        bool &finalize_topo, const std::string &xml_prefix="");
   void ReadNURBSMesh(std::istream &input, int &curved, int &read_gf,
                      bool spacing=false);
   void ReadInlineMesh(std::istream &input, bool generate_edges = false);
   void ReadGmshMesh(std::istream &input, int &curved, int &read_gf);

   /// Determine the mesh generator bitmask #meshgen, see MeshGenerator().
   /** Also, initializes #mesh_geoms. */
   void SetMeshGen();

   /// Return the length of the segment from node i to node j.
   real_t GetLength(int i, int j) const;

   void MarkForRefinement();
   void MarkTriMeshForRefinement();
   void GetEdgeOrdering(const DSTable &v_to_v, Array<int> &order);
   virtual void MarkTetMeshForRefinement(const DSTable &v_to_v);

   STable3D *GetFacesTable();
   STable3D *GetElementToFaceTable(int ret_ftbl = 0);

   /** Red refinement. Element with index i is refined. The default
       red refinement for now is Uniform. */
   void RedRefinement(int i, const DSTable &v_to_v,
                      int *edge1, int *edge2, int *middle)
   { UniformRefinement(i, v_to_v, edge1, edge2, middle); }

   /** Green refinement. Element with index i is refined. The default
       refinement for now is Bisection. */
   void GreenRefinement(int i, const DSTable &v_to_v,
                        int *edge1, int *edge2, int *middle)
   { Bisection(i, v_to_v, edge1, edge2, middle); }

   /// Bisect a triangle: element with index @a i is bisected.
   void Bisection(int i, const DSTable &, int *, int *, int *);

   /// Bisect a tetrahedron: element with index @a i is bisected.
   void Bisection(int i, HashTable<Hashed2> &);

   /// Bisect a boundary triangle: boundary element with index @a i is bisected.
   void BdrBisection(int i, const HashTable<Hashed2> &);

   /** Uniform Refinement. Element with index i is refined uniformly. */
   void UniformRefinement(int i, const DSTable &, int *, int *, int *);

   /// Returns the orientation of "test" relative to "base"
   static int GetTriOrientation(const int *base, const int *test);

   /// Returns the orientation of "base" relative to "test"
   /// In other words: GetTriOrientation(test, base) should equal
   /// InvertTriOrientation(GetTriOrientation(base, test))
   static int InvertTriOrientation(int ori);

   /// Returns the orientation of "c" relative to "a" by composing
   /// previously computed orientations relative to an intermediate
   /// set "b".
   static int ComposeTriOrientations(int ori_a_b, int ori_b_c);

   /// Returns the orientation of "test" relative to "base"
   static int GetQuadOrientation(const int *base, const int *test);

   /// Returns the orientation of "base" relative to "test"
   /// In other words: GetQuadOrientation(test, base) should equal
   /// InvertQuadOrientation(GetQuadOrientation(base, test))
   static int InvertQuadOrientation(int ori);

   /// Returns the orientation of "c" relative to "a" by composing
   /// previously computed orientations relative to an intermediate
   /// set "b".
   static int ComposeQuadOrientations(int ori_a_b, int ori_b_c);

   /// Returns the orientation of "test" relative to "base"
   static int GetTetOrientation(const int *base, const int *test);

   static void GetElementArrayEdgeTable(const Array<Element*> &elem_array,
                                        const DSTable &v_to_v,
                                        Table &el_to_edge);

   /** Return element to edge table and the indices for the boundary edges.
       The entries in the table are ordered according to the order of the
       nodes in the elements. For example, if T is the element to edge table
       T(i, 0) gives the index of edge in element i that connects vertex 0
       to vertex 1, etc. Returns the number of the edges. */
   int GetElementToEdgeTable(Table &);

   /// Used in GenerateFaces()
   void AddPointFaceElement(int lf, int gf, int el);

   void AddSegmentFaceElement (int lf, int gf, int el, int v0, int v1);

   void AddTriangleFaceElement (int lf, int gf, int el,
                                int v0, int v1, int v2);

   void AddQuadFaceElement (int lf, int gf, int el,
                            int v0, int v1, int v2, int v3);
   /** For a serial Mesh, return true if the face is interior. For a parallel
       ParMesh return true if the face is interior or shared. In parallel, this
       method only works if the face neighbor data is exchanged. */
   bool FaceIsTrueInterior(int FaceNo) const
   {
      return FaceIsInterior(FaceNo) || (faces_info[FaceNo].Elem2Inf >= 0);
   }

   void FreeElement(Element *E);

   void GenerateFaces();
   void GenerateNCFaceInfo();

   /// Begin construction of a mesh
   void InitMesh(int Dim_, int spaceDim_, int NVert, int NElem, int NBdrElem);

   // Used in the methods FinalizeXXXMesh() and FinalizeTopology()
   void FinalizeCheck();

   void Loader(std::istream &input, int generate_edges = 0,
               std::string parse_tag = "");

   /** If NURBS mesh, write NURBS format. If NCMesh, write mfem v1.1 format.
       If section_delimiter is empty, write mfem v1.0 format. Otherwise, write
       mfem v1.2 format with the given section_delimiter at the end.
       If @a comments is non-empty, it will be printed after the first line of
       the file, and each line should begin with '#'. */
   void Printer(std::ostream &os = mfem::out,
                std::string section_delimiter = "",
                const std::string &comments = "") const;

   /// @brief Creates a mesh for the parallelepiped [0,sx]x[0,sy]x[0,sz],
   /// divided into nx*ny*nz hexahedra if @a type = HEXAHEDRON or into
   /// 6*nx*ny*nz tetrahedrons if @a type = TETRAHEDRON.
   ///
   /// The parameter @a sfc_ordering controls how the elements
   /// (when @a type = HEXAHEDRON) are ordered: true - use space-filling curve
   /// ordering, or false - use lexicographic ordering.
   void Make3D(int nx, int ny, int nz, Element::Type type,
               real_t sx, real_t sy, real_t sz, bool sfc_ordering);

   /// @brief Creates a mesh for the parallelepiped [0,sx]x[0,sy]x[0,sz],
   /// divided into nx*ny*nz*24 tetrahedrons.
   ///
   /// The mesh is generated by taking nx*ny*nz hexahedra and splitting each
   /// hexahedron into 24 tetrahedrons. Each face of the hexahedron is split
   /// into 4 triangles (face edges are connected to a face-centered point),
   /// and the triangles are connected to a hex-centered point.
   void Make3D24TetsFromHex(int nx, int ny, int nz,
                            real_t sx, real_t sy, real_t sz);

   /// @brief Creates mesh for the rectangle [0,sx]x[0,sy], divided into nx*ny*4
   /// triangles.
   ///
   /// The mesh is generated by taking nx*ny quadrilaterals and splitting each
   /// quadrilateral into 4 triangles by connecting the vertices to a
   /// quad-centered point.
   void Make2D4TrisFromQuad(int nx, int ny, real_t sx, real_t sy);

   /// @brief Creates mesh for the rectangle [0,sx]x[0,sy], divided into nx*ny*5
   /// quadrilaterals.
   ///
   /// The mesh is generated by taking nx*ny quadrilaterals and splitting
   /// each quadrilateral into 5 quadrilaterals. Each quadrilateral is projected
   /// inwards and connected to the original quadrilateral.
   void Make2D5QuadsFromQuad(int nx, int ny, real_t sx, real_t sy);

   /// @brief Creates mesh for the rectangle [0,sx]x[0,sy], divided into nx*ny
   /// quadrilaterals if @a type = QUADRILATERAL or into 2*nx*ny triangles if
   /// @a type = TRIANGLE.
   ///
   /// If generate_edges = 0 (default) edges are not generated, if 1 edges are
   /// generated. The parameter @a sfc_ordering controls how the elements (when
   /// @a type = QUADRILATERAL) are ordered: true - use space-filling curve
   /// ordering, or false - use lexicographic ordering.
   void Make2D(int nx, int ny, Element::Type type, real_t sx, real_t sy,
               bool generate_edges, bool sfc_ordering);

   /// @a brief Creates a 1D mesh for the interval [0,sx] divided into n equal
   /// intervals.
   void Make1D(int n, real_t sx = 1.0);

   // used in GetElementData() and GetBdrElementData()
   void GetElementData(const Array<Element*> &elem_array, int geom,
                       Array<int> &elem_vtx, Array<int> &attr) const;

   // Internal helper used in MakeSimplicial (and ParMesh::MakeSimplicial).
   void MakeSimplicial_(const Mesh &orig_mesh, int *vglobal);

public:

   /// @anchor mfem_Mesh_ctors
   /// @name Standard Mesh constructors and related methods
   ///
   /// These constructors and assignment operators accept mesh information in
   /// a variety of common forms. For more specialized constructors see
   /// @ref mfem_Mesh_named_ctors "Named mesh constructors".
   /// @{
   Mesh() : attribute_sets(attributes), bdr_attribute_sets(bdr_attributes)
   { SetEmpty(); }

   /** Copy constructor. Performs a deep copy of (almost) all data, so that the
       source mesh can be modified (e.g. deleted, refined) without affecting the
       new mesh. If 'copy_nodes' is false, use a shallow (pointer) copy for the
       nodes, if present. */
   explicit Mesh(const Mesh &mesh, bool copy_nodes = true);

   /// Move constructor, useful for using a Mesh as a function return value.
   Mesh(Mesh &&mesh);

   /// Move assignment operator.
   Mesh& operator=(Mesh &&mesh);

   /// Explicitly delete the copy assignment operator.
   Mesh& operator=(const Mesh &mesh) = delete;

   /// Construct a Mesh from the given primary data.
   /** The array @a vertices is used as external data, i.e. the Mesh does not
       copy the data and will not delete the pointer.

       The data from the other arrays is copied into the internal Mesh data
       structures.

       This method calls the method FinalizeTopology(). The method Finalize()
       may be called after this constructor and after optionally setting the
       Mesh nodes. */
   Mesh(real_t *vertices, int num_vertices,
        int *element_indices, Geometry::Type element_type,
        int *element_attributes, int num_elements,
        int *boundary_indices, Geometry::Type boundary_type,
        int *boundary_attributes, int num_boundary_elements,
        int dimension, int space_dimension = -1);

   /** @anchor mfem_Mesh_init_ctor
       @brief _Init_ constructor: begin the construction of a Mesh object.

       Construct a shell of a mesh object allocating space to store pointers to
       the vertices, elements, and boundary elements. The vertices and elements
       themselves can later be added using methods from the
       @ref mfem_Mesh_construction "Mesh construction" group. */
   Mesh(int Dim_, int NVert, int NElem, int NBdrElem = 0, int spaceDim_ = -1)
      : attribute_sets(attributes), bdr_attribute_sets(bdr_attributes)
   {
      if (spaceDim_ == -1) { spaceDim_ = Dim_; }
      InitMesh(Dim_, spaceDim_, NVert, NElem, NBdrElem);
   }

   /** Creates mesh by reading a file in MFEM, Netgen, or VTK format. If
       generate_edges = 0 (default) edges are not generated, if 1 edges are
       generated. See also @a Mesh::LoadFromFile. See @a Mesh::Finalize for the
       meaning of @a refine. */
   explicit Mesh(const std::string &filename, int generate_edges = 0,
                 int refine = 1, bool fix_orientation = true);

   /** Creates mesh by reading data stream in MFEM, Netgen, or VTK format. If
       generate_edges = 0 (default) edges are not generated, if 1 edges are
       generated. */
   explicit Mesh(std::istream &input, int generate_edges = 0, int refine = 1,
                 bool fix_orientation = true);

   /// Create a disjoint mesh from the given mesh array
   ///
   /// @note Data is copied from the meshes in @a mesh_array.
   Mesh(Mesh *mesh_array[], int num_pieces);

   /** This is similar to the mesh constructor with the same arguments, but here
       the current mesh is destroyed and another one created based on the data
       stream again given in MFEM, Netgen, or VTK format. If generate_edges = 0
       (default) edges are not generated, if 1 edges are generated. */
   /// \see mfem::ifgzstream() for on-the-fly decompression of compressed ascii
   /// inputs.
   virtual void Load(std::istream &input, int generate_edges = 0,
                     int refine = 1, bool fix_orientation = true)
   {
      Loader(input, generate_edges);
      Finalize(refine, fix_orientation);
   }

   /// Swaps internal data with another mesh. By default, non-geometry members
   /// like 'ncmesh' and 'NURBSExt' are only swapped when 'non_geometry' is set.
   void Swap(Mesh& other, bool non_geometry);

   /// Clear the contents of the Mesh.
   void Clear() { Destroy(); SetEmpty(); }

   /// Destroys Mesh.
   virtual ~Mesh() { DestroyPointers(); }

   /// @}

   /** @anchor mfem_Mesh_named_ctors @name Named mesh constructors.

        Each of these constructors uses the move constructor, and can be used as
        the right-hand side of an assignment when creating new meshes. For more
        general mesh constructors see
        @ref mfem_Mesh_ctors "Standard mesh constructors".*/
   ///@{

   /** Creates mesh by reading a file in MFEM, Netgen, or VTK format. If
       generate_edges = 0 (default) edges are not generated, if 1 edges are
       generated.

       @note @a filename is not cached by the Mesh object and can be
       safely deleted following this function call.
   */
   static Mesh LoadFromFile(const std::string &filename,
                            int generate_edges = 0, int refine = 1,
                            bool fix_orientation = true);

   /// Creates 1D mesh, divided into n equal intervals.
   static Mesh MakeCartesian1D(int n, real_t sx = 1.0);

   /// @brief Creates mesh for the rectangle [0,sx]x[0,sy], divided into nx*ny
   /// quadrilaterals if @a type = QUADRILATERAL or into 2*nx*ny triangles if
   /// @a type = TRIANGLE.
   ///
   /// If generate_edges = 0 (default) edges are not generated, if 1 edges are
   /// generated. The parameter @a sfc_ordering controls how the elements (when
   /// @a type = QUADRILATERAL) are ordered: true - use space-filling curve
   /// ordering, or false - use lexicographic ordering.
   static Mesh MakeCartesian2D(
      int nx, int ny, Element::Type type, bool generate_edges = false,
      real_t sx = 1.0, real_t sy = 1.0, bool sfc_ordering = true);

   /// @brief Creates a mesh for the parallelepiped [0,sx]x[0,sy]x[0,sz],
   /// divided into nx*ny*nz hexahedra if @a type = HEXAHEDRON or into
   /// 6*nx*ny*nz tetrahedrons if @a type = TETRAHEDRON.
   ///
   /// The parameter @a sfc_ordering controls how the elements
   /// (when @a type = HEXAHEDRON) are ordered: true - use space-filling curve
   /// ordering, or false - use lexicographic ordering.
   static Mesh MakeCartesian3D(
      int nx, int ny, int nz, Element::Type type,
      real_t sx = 1.0, real_t sy = 1.0, real_t sz = 1.0,
      bool sfc_ordering = true);

   /// @brief Creates a mesh for the parallelepiped [0,sx]x[0,sy]x[0,sz],
   /// divided into nx*ny*nz*24 tetrahedrons.
   ///
   /// The mesh is generated by taking nx*ny*nz hexahedra and splitting each
   /// hexahedron into 24 tetrahedrons. Each face of the hexahedron is split
   /// into 4 triangles (face edges are connected to a face-centered point),
   /// and the triangles are connected to a hex-centered point.
   static Mesh MakeCartesian3DWith24TetsPerHex(int nx, int ny, int nz,
                                               real_t sx = 1.0, real_t sy = 1.0,
                                               real_t sz = 1.0);

   /// @brief Creates mesh for the rectangle [0,sx]x[0,sy], divided into nx*ny*4
   /// triangles.
   ///
   /// The mesh is generated by taking nx*ny quadrilaterals and splitting each
   /// quadrilateral into 4 triangles by connecting the vertices to a
   /// quad-centered point.
   static Mesh MakeCartesian2DWith4TrisPerQuad(int nx, int ny, real_t sx = 1.0,
                                               real_t sy = 1.0);

   /// @brief Creates mesh for the rectangle [0,sx]x[0,sy], divided into nx*ny*5
   /// quadrilaterals.
   ///
   /// The mesh is generated by taking nx*ny quadrilaterals and splitting
   /// each quadrilateral into 5 quadrilaterals. Each quadrilateral is projected
   /// inwards and connected to the original quadrilateral.
   static Mesh MakeCartesian2DWith5QuadsPerQuad(int nx, int ny, real_t sx = 1.0,
                                                real_t sy = 1.0);


   /// Create a refined (by any factor) version of @a orig_mesh.
   /** @param[in] orig_mesh  The starting coarse mesh.
       @param[in] ref_factor The refinement factor, an integer > 1.
       @param[in] ref_type   Specify the positions of the new vertices. The
                             options are BasisType::ClosedUniform or
                             BasisType::GaussLobatto.

       The refinement data which can be accessed with GetRefinementTransforms()
       is set to reflect the performed refinements.

       @note The constructed Mesh is straight-sided. */
   static Mesh MakeRefined(Mesh &orig_mesh, int ref_factor, int ref_type);

   /// Create a refined mesh, where each element of the original mesh may be
   /// refined by a different factor.
   /** @param[in] orig_mesh   The starting coarse mesh.
       @param[in] ref_factors An array of integers whose size is the number of
                              elements of @a orig_mesh. The @a ith element of
                              @a orig_mesh is refined by refinement factor
                              @a ref_factors[i].
       @param[in] ref_type    Specify the positions of the new vertices. The
                              options are BasisType::ClosedUniform or
                              BasisType::GaussLobatto.

       The refinement data which can be accessed with GetRefinementTransforms()
       is set to reflect the performed refinements.

       @note The constructed Mesh is straight-sided. */
   /// refined @a ref_factors[i] times in each dimension.
   static Mesh MakeRefined(Mesh &orig_mesh, const Array<int> &ref_factors,
                           int ref_type);

   /** Create a mesh by splitting each element of @a orig_mesh into simplices.
       Quadrilaterals are split into two triangles, prisms are split into
       3 tetrahedra, and hexahedra are split into either 5 or 6 tetrahedra
       depending on the configuration.
       @warning The curvature of the original mesh is not carried over to the
       new mesh. Periodic meshes are not supported. */
   static Mesh MakeSimplicial(const Mesh &orig_mesh);

   /// Create a periodic mesh by identifying vertices of @a orig_mesh.
   /** Each vertex @a i will be mapped to vertex @a v2v[i], such that all
       vertices that are coincident under the periodic mapping get mapped to
       the same index. The mapping @a v2v can be generated from translation
       vectors using Mesh::CreatePeriodicVertexMapping.
       @note MFEM requires that each edge of the resulting mesh be uniquely
       identifiable by a pair of distinct vertices. As a consequence, periodic
       boundaries must be separated by at least two interior vertices.
       @note The resulting mesh uses a discontinuous nodal function, see
       SetCurvature() for further details. */
   static Mesh MakePeriodic(const Mesh &orig_mesh, const std::vector<int> &v2v);

   ///@}

   /** @anchor mfem_Mesh_construction
       @name Methods for piecewise Mesh construction.

       These methods are intended to be used with the @ref mfem_Mesh_init_ctor
       "init constructor". */
   ///@{

   /// @note The returned object should be deleted by the caller.
   Element *NewElement(int geom);

   int AddVertex(real_t x, real_t y = 0.0, real_t z = 0.0);
   int AddVertex(const real_t *coords);
   int AddVertex(const Vector &coords);
   /// Mark vertex @a i as nonconforming, with parent vertices @a p1 and @a p2.
   void AddVertexParents(int i, int p1, int p2);
   /// Adds a vertex at the mean center of the @a nverts vertex indices given
   /// by @a vi.
   int AddVertexAtMeanCenter(const int *vi, const int nverts, int dim = 3);

   /// Adds a segment to the mesh given by 2 vertices @a v1 and @a v2.
   int AddSegment(int v1, int v2, int attr = 1);
   /// Adds a segment to the mesh given by 2 vertices @a vi.
   int AddSegment(const int *vi, int attr = 1);

   /// Adds a triangle to the mesh given by 3 vertices @a v1 through @a v3.
   int AddTriangle(int v1, int v2, int v3, int attr = 1);
   /// Adds a triangle to the mesh given by 3 vertices @a vi.
   int AddTriangle(const int *vi, int attr = 1);
   /// Adds a triangle to the mesh given by 3 vertices @a vi.
   int AddTri(const int *vi, int attr = 1) { return AddTriangle(vi, attr); }

   /// Adds a quadrilateral to the mesh given by 4 vertices @a v1 through @a v4.
   int AddQuad(int v1, int v2, int v3, int v4, int attr = 1);
   /// Adds a quadrilateral to the mesh given by 4 vertices @a vi.
   int AddQuad(const int *vi, int attr = 1);

   /// Adds a tetrahedron to the mesh given by 4 vertices @a v1 through @a v4.
   int AddTet(int v1, int v2, int v3, int v4, int attr = 1);
   /// Adds a tetrahedron to the mesh given by 4 vertices @a vi.
   int AddTet(const int *vi, int attr = 1);

   /// Adds a wedge to the mesh given by 6 vertices @a v1 through @a v6.
   int AddWedge(int v1, int v2, int v3, int v4, int v5, int v6, int attr = 1);
   /// Adds a wedge to the mesh given by 6 vertices @a vi.
   int AddWedge(const int *vi, int attr = 1);

   /// Adds a pyramid to the mesh given by 5 vertices @a v1 through @a v5.
   int AddPyramid(int v1, int v2, int v3, int v4, int v5, int attr = 1);
   /// Adds a pyramid to the mesh given by 5 vertices @a vi.
   int AddPyramid(const int *vi, int attr = 1);

   /// Adds a hexahedron to the mesh given by 8 vertices @a v1 through @a v8.
   int AddHex(int v1, int v2, int v3, int v4, int v5, int v6, int v7, int v8,
              int attr = 1);
   /// Adds a hexahedron to the mesh given by 8 vertices @a vi.
   int AddHex(const int *vi, int attr = 1);
   /// @brief Adds 6 tetrahedrons to the mesh by splitting a hexahedron given by
   /// 8 vertices @a vi.
   void AddHexAsTets(const int *vi, int attr = 1);
   /// @brief Adds 2 wedges to the mesh by splitting a hexahedron given by
   /// 8 vertices @a vi.
   void AddHexAsWedges(const int *vi, int attr = 1);
   /// @brief Adds 6 pyramids to the mesh by splitting a hexahedron given by
   /// 8 vertices @a vi.
   void AddHexAsPyramids(const int *vi, int attr = 1);

   /// @brief Adds 24 tetrahedrons to the mesh by splitting a hexahedron.
   ///
   /// @a vi are the 8 vertices of the hexahedron, @a hex_face_verts has the
   /// map from the 4 vertices of each face of the hexahedron to the index
   /// of the point created at the center of the face, and @a attr is the
   /// attribute of the new elements. See @a Make3D24TetsFromHex for usage.
   void AddHexAs24TetsWithPoints(int *vi,
                                 std::map<std::array<int, 4>, int>
                                 &hex_face_verts,
                                 int attr = 1);

   /// @brief Adds 4 triangles to the mesh by splitting a quadrilateral given by
   /// 4 vertices @a vi.
   ///
   /// @a attr is the attribute of the new elements. See @a Make2D4TrisFromQuad
   /// for usage.
   void AddQuadAs4TrisWithPoints(int *vi, int attr = 1);

   /// @brief Adds 5 quadrilaterals to the mesh by splitting a quadrilateral
   /// given by 4 vertices @a vi.
   ///
   /// @a attr is the attribute of the new elements. See @a Make2D5QuadsFromQuad
   /// for usage.
   void AddQuadAs5QuadsWithPoints(int *vi, int attr = 1);

   /// The parameter @a elem should be allocated using the NewElement() method
   /// @note Ownership of @a elem will pass to the Mesh object
   int AddElement(Element *elem);
   /// The parameter @a elem should be allocated using the NewElement() method
   /// @note Ownership of @a elem will pass to the Mesh object
   int AddBdrElement(Element *elem);

   /**
    * @brief Add an array of boundary elements to the mesh, along with map from
    * the elements to their faces
    * @param[in] bdr_elems The set of boundary element pointers, ownership of
    * the pointers will be transferred to the Mesh object
    * @param[in] be_to_face The map from the boundary element index to the face
    * index
    */
   void AddBdrElements(Array<Element *> &bdr_elems,
                       const Array<int> &be_to_face);

   int AddBdrSegment(int v1, int v2, int attr = 1);
   int AddBdrSegment(const int *vi, int attr = 1);

   int AddBdrTriangle(int v1, int v2, int v3, int attr = 1);
   int AddBdrTriangle(const int *vi, int attr = 1);

   int AddBdrQuad(int v1, int v2, int v3, int v4, int attr = 1);
   int AddBdrQuad(const int *vi, int attr = 1);
   void AddBdrQuadAsTriangles(const int *vi, int attr = 1);

   int AddBdrPoint(int v, int attr = 1);

   virtual void GenerateBoundaryElements();
   /// Finalize the construction of a triangular Mesh.
   void FinalizeTriMesh(int generate_edges = 0, int refine = 0,
                        bool fix_orientation = true);
   /// Finalize the construction of a quadrilateral Mesh.
   void FinalizeQuadMesh(int generate_edges = 0, int refine = 0,
                         bool fix_orientation = true);
   /// Finalize the construction of a tetrahedral Mesh.
   void FinalizeTetMesh(int generate_edges = 0, int refine = 0,
                        bool fix_orientation = true);
   /// Finalize the construction of a wedge Mesh.
   void FinalizeWedgeMesh(int generate_edges = 0, int refine = 0,
                          bool fix_orientation = true);
   /// Finalize the construction of a hexahedral Mesh.
   void FinalizeHexMesh(int generate_edges = 0, int refine = 0,
                        bool fix_orientation = true);
   /// Finalize the construction of any type of Mesh.
   /** This method calls FinalizeTopology() and Finalize(). */
   void FinalizeMesh(int refine = 0, bool fix_orientation = true);

   ///@}

   /// @name Mesh consistency methods
   /// @{

   /** @brief Finalize the construction of the secondary topology (connectivity)
       data of a Mesh. */
   /** This method does not require any actual coordinate data (either vertex
       coordinates for linear meshes or node coordinates for meshes with nodes)
       to be available. However, the data generated by this method is generally
       required by the FiniteElementSpace class.

       After calling this method, setting the Mesh vertices or nodes, it may be
       appropriate to call the method Finalize(). */
   void FinalizeTopology(bool generate_bdr = true);

   /// Finalize the construction of a general Mesh.
   /** This method will:
       - check and optionally fix the orientation of regular elements
       - check and fix the orientation of boundary elements
       - assume that #vertices are defined, if #Nodes == NULL
       - assume that #Nodes are defined, if #Nodes != NULL.
       @param[in] refine  If true, prepare the Mesh for conforming refinement of
                          triangular or tetrahedral meshes.
       @param[in] fix_orientation
                          If true, fix the orientation of inverted mesh elements
                          by permuting their vertices.

       Before calling this method, call FinalizeTopology() and ensure that the
       Mesh vertices or nodes are set. */
   virtual void Finalize(bool refine = false, bool fix_orientation = false);

   /// @brief Determine the sets of unique attribute values in domain and
   /// boundary elements.
   ///
   /// Separately scan the domain and boundary elements to generate unique,
   /// sorted sets of the element attribute values present in the mesh and
   /// store these in the Mesh::attributes and Mesh::bdr_attributes arrays.
   virtual void SetAttributes();

   /// Remove unused vertices and rebuild mesh connectivity.
   void RemoveUnusedVertices();

   /** Remove boundary elements that lie in the interior of the mesh, i.e. that
       have two adjacent faces in 3D, or edges in 2D. */
   void RemoveInternalBoundaries();

   /**
    * @brief Clear the boundary element to edge map.
    */
   void DeleteBoundaryElementToEdge()
   {
      delete bel_to_edge;
      bel_to_edge = nullptr;
   }

   /// @}

   /// @name Element ordering methods
   /// @{

   /** This is our integration with the Gecko library. The method finds an
       element ordering that will increase memory coherency by putting elements
       that are in physical proximity closer in memory. It can also be used to
       obtain a space-filling curve ordering for ParNCMesh partitioning.
       @param[out] ordering Output element ordering.
       @param iterations Total number of V cycles. The ordering may improve with
       more iterations. The best iteration is returned at the end.
       @param window Initial window size. This determines the number of
       permutations tested at each multigrid level and strongly influences the
       quality of the result, but the cost of increasing 'window' is exponential.
       @param period The window size is incremented every 'period' iterations.
       @param seed Seed for initial random ordering (0 = skip random reorder).
       @param verbose Print the progress of the optimization to mfem::out.
       @param time_limit Optional time limit for the optimization, in seconds.
       When reached, ordering from the best iteration so far is returned
       (0 = no limit).
       @return The final edge product cost of the ordering. The function may be
       called in an external loop with different seeds, and the best ordering can
       then be retained. */
   real_t GetGeckoElementOrdering(Array<int> &ordering,
                                  int iterations = 4, int window = 4,
                                  int period = 2, int seed = 0,
                                  bool verbose = false, real_t time_limit = 0);

   /** Return an ordering of the elements that approximately follows the Hilbert
       curve. The method performs a spatial (Hilbert) sort on the centers of all
       elements and returns the resulting sequence, which can then be passed to
       ReorderElements. This is a cheap alternative to GetGeckoElementOrdering.*/
   void GetHilbertElementOrdering(Array<int> &ordering);

   /** Rebuilds the mesh with a different order of elements. For each element i,
       the array ordering[i] contains its desired new index. Note that the method
       reorders vertices, edges and faces along with the elements. */
   void ReorderElements(const Array<int> &ordering, bool reorder_vertices = true);

   /// @}

   /// @name Information about the mesh as a whole
   /// @{

   /// @brief Dimension of the reference space used within the elements
   int Dimension() const { return Dim; }

   /// @brief Dimension of the physical space containing the mesh
   int SpaceDimension() const { return spaceDim; }

   /// Equals 1 + num_holes - num_loops
   inline int EulerNumber() const
   { return NumOfVertices - NumOfEdges + NumOfFaces - NumOfElements; }
   /// Equals 1 - num_holes
   inline int EulerNumber2D() const
   { return NumOfVertices - NumOfEdges + NumOfElements; }

   /** @brief Get the mesh generator/type.

       The purpose of this is to be able to quickly tell what type of elements
       one has in the mesh. Examination of this bitmask along with knowledge
       of the mesh dimension can be used to identify which element types are
       present.

       @return A bitmask:
       - bit 0 - simplices are present in the mesh (triangles, tets),
       - bit 1 - tensor product elements are present in the mesh (quads, hexes),
       - bit 2 - the mesh has wedge elements.
       - bit 3 - the mesh has pyramid elements.

       In parallel, the result takes into account elements on all processors.
   */
   inline int MeshGenerator() const { return meshgen; }

   /// Checks if the mesh has boundary elements
   virtual bool HasBoundaryElements() const { return (NumOfBdrElements > 0); }

   /** @brief Return true iff the given @a geom is encountered in the mesh.
       Geometries of dimensions lower than Dimension() are counted as well. */
   bool HasGeometry(Geometry::Type geom) const
   { return mesh_geoms & (1 << geom); }

   /** @brief Return the number of geometries of the given dimension present in
       the mesh. */
   /** For a parallel mesh only the local geometries are counted. */
   int GetNumGeometries(int dim) const;

   /// Return all element geometries of the given dimension present in the mesh.
   /** For a parallel mesh only the local geometries are returned.

       The returned geometries are sorted. */
   void GetGeometries(int dim, Array<Geometry::Type> &el_geoms) const;

   void GetCharacteristics(real_t &h_min, real_t &h_max,
                           real_t &kappa_min, real_t &kappa_max,
                           Vector *Vh = NULL, Vector *Vk = NULL);

   /// @}

   /// @name Information concerning numbers of mesh entities
   /// @{

   /** @brief Returns number of vertices.  Vertices are only at the corners of
       elements, where you would expect them in the lowest-order mesh. */
   inline int GetNV() const { return NumOfVertices; }

   /// Returns number of elements.
   inline int GetNE() const { return NumOfElements; }

   /// Returns number of boundary elements.
   inline int GetNBE() const { return NumOfBdrElements; }

   /// Return the number of edges.
   inline int GetNEdges() const { return NumOfEdges; }

   /// Return the number of faces in a 3D mesh.
   inline int GetNFaces() const { return NumOfFaces; }

   /// Return the number of faces (3D), edges (2D) or vertices (1D).
   int GetNumFaces() const;

   /** @brief Return the number of faces (3D), edges (2D) or vertices (1D)
       including ghost faces. */
   int GetNumFacesWithGhost() const;

   /** @brief Returns the number of faces according to the requested type, does
       not count master nonconforming faces.

       If type==Boundary returns only the number of true boundary faces
       contrary to GetNBE() that returns all "boundary" elements which may
       include actual interior faces.
       Similarly, if type==Interior, only the true interior faces are counted
       excluding all master nonconforming faces. */
   virtual int GetNFbyType(FaceType type) const;

   /// Return the total (global) number of elements.
   long long GetGlobalNE() const { return NumOfElements; }

   /// @}

   /// @name Access to individual mesh entities
   /// @{

   /// @brief Return pointer to vertex i's coordinates.
   /// @warning For high-order meshes (when #Nodes != NULL) vertices may not be
   /// updated and should not be used!
   const real_t *GetVertex(int i) const { return vertices[i](); }

   /// @brief Return pointer to vertex i's coordinates.
   ///
   /// @warning For high-order meshes (when Nodes != NULL) vertices may not
   /// being updated and should not be used!
   ///
   /// @note The pointer returned by this function can be used to
   /// alter vertex locations but the pointer itself should not be
   /// changed by the caller.
   real_t *GetVertex(int i) { return vertices[i](); }

   /// @brief Return pointer to the i'th element object
   ///
   /// The index @a i should be in the range [0, Mesh::GetNE())
   ///
   /// In parallel, @a i is the local element index which is in the
   /// same range mentioned above.
   const Element *GetElement(int i) const { return elements[i]; }

   /// @brief Return pointer to the i'th element object
   ///
   /// @note Provides read/write access to the i'th element object so
   /// that element attributes or connectivity can be adjusted. However,
   /// the Element object itself should not be deleted by the caller.
   Element *GetElement(int i) { return elements[i]; }

   /// @brief Return pointer to the i'th boundary element object
   ///
   /// The index @a i should be in the range [0, Mesh::GetNBE())
   ///
   /// In parallel, @a i is the local boundary element index which is
   /// in the same range mentioned above.
   const Element *GetBdrElement(int i) const { return boundary[i]; }

   /// @brief Return pointer to the i'th boundary element object
   ///
   /// @note Provides read/write access to the i'th boundary element object so
   /// that boundary attributes or connectivity can be adjusted. However,
   /// the Element object itself should not be deleted by the caller.
   Element *GetBdrElement(int i) { return boundary[i]; }

   /// @brief Return pointer to the i'th face element object
   ///
   /// The index @a i should be in the range [0, Mesh::GetNFaces())
   const Element *GetFace(int i) const { return faces[i]; }

   /// @}

   /// @name Access to groups of mesh entities
   /// @{

   const Element* const *GetElementsArray() const
   { return elements.GetData(); }

   void GetElementData(int geom, Array<int> &elem_vtx, Array<int> &attr) const
   { GetElementData(elements, geom, elem_vtx, attr); }

   void GetBdrElementData(int geom, Array<int> &bdr_elem_vtx,
                          Array<int> &bdr_attr) const
   { GetElementData(boundary, geom, bdr_elem_vtx, bdr_attr); }

   /// @}

   /// @name Access information concerning individual mesh entites
   /// @{

   /// Return the attribute of element i.
   int GetAttribute(int i) const { return elements[i]->GetAttribute(); }

   /// Set the attribute of element i.
   void SetAttribute(int i, int attr);

   /// Return the attribute of boundary element i.
   int GetBdrAttribute(int i) const { return boundary[i]->GetAttribute(); }

   /// Set the attribute of boundary element i.
   void SetBdrAttribute(int i, int attr) { boundary[i]->SetAttribute(attr); }

   /// Return the attribute of patch i, for a NURBS mesh.
   int GetPatchAttribute(int i) const;

   /// Set the attribute of patch i, for a NURBS mesh.
   void SetPatchAttribute(int i, int attr);

   /// Return the attribute of patch boundary element i, for a NURBS mesh.
   int GetPatchBdrAttribute(int i) const;

   /// Set the attribute of patch boundary element i, for a NURBS mesh.
   void SetPatchBdrAttribute(int i, int attr);

   /// Returns the type of element i.
   Element::Type GetElementType(int i) const;

   /// Returns the type of boundary element i.
   Element::Type GetBdrElementType(int i) const;

   Element::Type  GetFaceElementType(int Face) const;

   /// Return the Geometry::Type associated with face @a i.
   Geometry::Type GetFaceGeometry(int i) const;

   Geometry::Type GetElementGeometry(int i) const
   {
      return elements[i]->GetGeometryType();
   }

   Geometry::Type GetBdrElementGeometry(int i) const
   {
      return boundary[i]->GetGeometryType();
   }

   Geometry::Type GetElementBaseGeometry(int i) const
   { return GetElementGeometry(i); }

   Geometry::Type GetBdrElementBaseGeometry(int i) const
   { return GetBdrElementGeometry(i); }

   /// Return true if the given face is interior. @sa FaceIsTrueInterior().
   bool FaceIsInterior(int FaceNo) const
   {
      return (faces_info[FaceNo].Elem2No >= 0);
   }

   /** @brief Get the size of the i-th element relative to the perfect
       reference element. */
   real_t GetElementSize(int i, int type = 0);

   real_t GetElementSize(int i, const Vector &dir);

   real_t GetElementVolume(int i);

   void GetElementCenter(int i, Vector &center);

   /** Compute the Jacobian of the transformation from the perfect
       reference element at the given integration point (defaults to the
       center of the element if no integration point is specified) */
   void GetElementJacobian(int i, DenseMatrix &J,
                           const IntegrationPoint *ip = NULL);

   /// @}

   /// List of mesh geometries stored as Array<Geometry::Type>.
   class GeometryList : public Array<Geometry::Type>
   {
   protected:
      Geometry::Type geom_buf[Geometry::NumGeom];
   public:
      /// Construct a GeometryList of all element geometries in @a mesh.
      GeometryList(const Mesh &mesh)
         : Array<Geometry::Type>(geom_buf, Geometry::NumGeom)
      { mesh.GetGeometries(mesh.Dimension(), *this); }
      /** @brief Construct a GeometryList of all geometries of dimension @a dim
          in @a mesh. */
      GeometryList(const Mesh &mesh, int dim)
         : Array<Geometry::Type>(geom_buf, Geometry::NumGeom)
      { mesh.GetGeometries(dim, *this); }
   };

   /// @name Access connectivity for individual mesh entites
   /// @{

   /// Returns the indices of the vertices of element i.
   void GetElementVertices(int i, Array<int> &v) const
   { elements[i]->GetVertices(v); }

   /// Returns the indices of the vertices of boundary element i.
   void GetBdrElementVertices(int i, Array<int> &v) const
   { boundary[i]->GetVertices(v); }

   /// Return the indices and the orientations of all edges of element i.
   void GetElementEdges(int i, Array<int> &edges, Array<int> &cor) const;

   /// Return the indices and the orientations of all edges of bdr element i.
   void GetBdrElementEdges(int i, Array<int> &edges, Array<int> &cor) const;

   /** Return the indices and the orientations of all edges of face i.
       Works for both 2D (face=edge) and 3D faces. */
   void GetFaceEdges(int i, Array<int> &edges, Array<int> &o) const;

   /// Returns the indices of the vertices of face i.
   void GetFaceVertices(int i, Array<int> &vert) const
   {
      if (Dim == 1)
      {
         vert.SetSize(1); vert[0] = i;
      }
      else
      {
         faces[i]->GetVertices(vert);
      }
   }

   /// Returns the indices of the vertices of edge i.
   void GetEdgeVertices(int i, Array<int> &vert) const;

   /// Return the indices and the orientations of all faces of element i.
   void GetElementFaces(int i, Array<int> &faces, Array<int> &ori) const;

   /** @brief Returns the sorted, unique indices of elements sharing a face with
       element @a elem, including @a elem. */
   Array<int> FindFaceNeighbors(const int elem) const;

   /** Return the index and the orientation of the vertex of bdr element i. (1D)
       Return the index and the orientation of the edge of bdr element i. (2D)
       Return the index and the orientation of the face of bdr element i. (3D)

       In 2D, the returned edge orientation is 0 or 1, not +/-1 as returned by
       GetElementEdges/GetBdrElementEdges. */
   void GetBdrElementFace(int i, int *f, int *o) const;

   /** @brief For the given boundary element, bdr_el, return its adjacent
       element and its info, i.e. 64*local_bdr_index+bdr_orientation.

       The returned bdr_orientation is that of the boundary element relative to
       the respective face element.

       @sa GetBdrElementAdjacentElement2() */
   void GetBdrElementAdjacentElement(int bdr_el, int &el, int &info) const;

   /// @brief Return the local face (codimension-1) index for the given boundary
   /// element index.
   int GetBdrElementFaceIndex(int be_idx) const { return be_to_face[be_idx]; }

   /// @}

   /// @name Access connectivity data
   /// @{

   /// @note The returned Table should be deleted by the caller
   Table *GetVertexToElementTable();

   /// @note The returned Table should be deleted by the caller
   Table *GetVertexToBdrElementTable();

   /// Return the "face"-element Table. Here "face" refers to face (3D),
   /// edge (2D), or vertex (1D).
   ///
   /// @note The returned Table should be deleted by the caller.
   Table *GetFaceToElementTable() const;

   /// Returns the face-to-edge Table (3D)
   ///
   /// @note The returned object should NOT be deleted by the caller.
   Table *GetFaceEdgeTable() const;

   /// Returns the edge-to-vertex Table (3D)
   ///
   /// @note The returned object should NOT be deleted by the caller.
   Table *GetEdgeVertexTable() const;

   /** Return vertex to vertex table. The connections stored in the table
    are from smaller to bigger vertex index, i.e. if i<j and (i, j) is
    in the table, then (j, i) is not stored.

    @note This data is not stored internally as a Table. The Table passed as
    an argument is populated using the EdgeVertex Table (see GetEdgeVertexTable)
    if available or the element connectivity.
   */
   void GetVertexToVertexTable(DSTable &) const;

   const Table &ElementToElementTable();

   const Table &ElementToFaceTable() const;

   const Table &ElementToEdgeTable() const;

   Array<int> GetFaceToBdrElMap() const;

   ///@}

   /** This enumerated type describes the three main face topologies:
       - Boundary, for faces on the boundary of the computational domain,
       - Conforming, for conforming faces interior to the computational domain,
       - Nonconforming, for nonconforming faces interior to the computational
         domain. */
   enum class FaceTopology { Boundary,
                             Conforming,
                             Nonconforming,
                             NA
                           };

   /** This enumerated type describes the location of the two elements sharing a
       face, Local meaning that the element is local to the MPI rank, FaceNbr
       meaning that the element is distributed on a different MPI rank, this
       typically means that methods with FaceNbr should be used to access the
       relevant information, e.g., ParFiniteElementSpace::GetFaceNbrElementVDofs.
    */
   enum class ElementLocation { Local, FaceNbr, NA };

   /** This enumerated type describes the topological relation of an element to
       a face:
       - Coincident meaning that the element's face is topologically equal to
         the mesh face.
       - Superset meaning that the element's face is topologically coarser than
         the mesh face, i.e., the element's face contains the mesh face.
       - Subset meaning that the element's face is topologically finer than the
         mesh face, i.e., the element's face is contained in the mesh face.
       Superset and Subset are only relevant for nonconforming faces.
       Master nonconforming faces have a conforming element on one side, and a
       fine element on the other side. Slave nonconforming faces have a
       conforming element on one side, and a coarse element on the other side.
    */
   enum class ElementConformity { Coincident, Superset, Subset, NA };

   /** This enumerated type describes the corresponding FaceInfo internal
       representation (encoded cases), c.f. FaceInfo's documentation:
       Classification of a local (non-ghost) face based on its FaceInfo:
         - Elem2No >= 0 --> local interior face; can be either:
            - NCFace == -1 --> LocalConforming,
            - NCFace >= 0 --> LocalSlaveNonconforming,
         - Elem2No < 0 --> local "boundary" face; can be one of:
            - NCFace == -1 --> conforming face; can be either:
               - Elem2Inf < 0 --> Boundary,
               - Elem2Inf >= 0 --> SharedConforming,
            - NCFace >= 0 --> nonconforming face; can be one of:
               - Elem2Inf < 0 --> MasterNonconforming (shared or not shared),
               - Elem2Inf >= 0 --> SharedSlaveNonconforming.
       Classification of a ghost (non-local) face based on its FaceInfo:
         - Elem1No == -1 --> GhostMaster (includes other unused ghost faces),
         - Elem1No >= 0 --> GhostSlave.
    */
   enum class FaceInfoTag { Boundary,
                            LocalConforming,
                            LocalSlaveNonconforming,
                            SharedConforming,
                            SharedSlaveNonconforming,
                            MasterNonconforming,
                            GhostSlave,
                            GhostMaster
                          };

   /** @brief This structure is used as a human readable output format that
       deciphers the information contained in Mesh::FaceInfo when using the
       Mesh::GetFaceInformation() method.

       The element indices in this structure don't need further processing,
       contrary to the ones obtained through Mesh::GetFacesElements and can
       directly be used, e.g., Elem1 and Elem2 indices.
       Likewise the orientations for Elem1 and Elem2 already take into account
       special cases and can be used as is.
   */
   struct FaceInformation
   {
      FaceTopology topology;

      struct
      {
         ElementLocation location;
         ElementConformity conformity;
         int index;
         int local_face_id;
         int orientation;
      } element[2];

      FaceInfoTag tag;
      int ncface;
      const DenseMatrix* point_matrix;

      /** @brief Return true if the face is a local interior face which is NOT
          a master nonconforming face. */
      bool IsLocal() const
      {
         return element[1].location == Mesh::ElementLocation::Local;
      }

      /** @brief Return true if the face is a shared interior face which is NOT
          a master nonconforming face. */
      bool IsShared() const
      {
         return element[1].location == Mesh::ElementLocation::FaceNbr;
      }

      /** @brief return true if the face is an interior face to the computation
          domain, either a local or shared interior face (not a boundary face)
          which is NOT a master nonconforming face.
       */
      bool IsInterior() const
      {
         return topology == FaceTopology::Conforming ||
                topology == FaceTopology::Nonconforming;
      }

      /** @brief Return true if the face is a boundary face. */
      bool IsBoundary() const
      {
         return topology == FaceTopology::Boundary;
      }

      /// @brief Return true if the face is of the same type as @a type.
      bool IsOfFaceType(FaceType type) const
      {
         switch (type)
         {
            case FaceType::Interior:
               return IsInterior();
            case FaceType::Boundary:
               return IsBoundary();
            default:
               return false;
         }
      }

      /// @brief Return true if the face is a conforming face.
      bool IsConforming() const
      {
         return topology == FaceTopology::Conforming;
      }

      /// @brief Return true if the face is a nonconforming fine face.
      bool IsNonconformingFine() const
      {
         return topology == FaceTopology::Nonconforming &&
                (element[0].conformity == ElementConformity::Superset ||
                 element[1].conformity == ElementConformity::Superset);
      }

      /// @brief Return true if the face is a nonconforming coarse face.
      /** Note that ghost nonconforming master faces cannot be clearly
          identified as such with the currently available information, so this
          method will return false for such faces. */
      bool IsNonconformingCoarse() const
      {
         return topology == FaceTopology::Nonconforming &&
                element[1].conformity == ElementConformity::Subset;
      }

      /// @brief cast operator from FaceInformation to FaceInfo.
      operator Mesh::FaceInfo() const;
   };

   /// Given a "face info int", return the face orientation. @sa FaceInfo.
   static int DecodeFaceInfoOrientation(int info) { return info%64; }

   /// Given a "face info int", return the local face index. @sa FaceInfo.
   static int DecodeFaceInfoLocalIndex(int info) { return info/64; }

   /// @brief Given @a local_face_index and @a orientation, return the
   /// corresponding encoded "face info int". @sa FaceInfo.
   static int EncodeFaceInfo(int local_face_index, int orientation)
   { return orientation + local_face_index*64; }

   /// @name More advanced entity information access methods
   /// @{

   /* Return point matrix of element i of dimension Dim X #v, where for every
      vertex we give its coordinates in space of dimension Dim. */
   void GetPointMatrix(int i, DenseMatrix &pointmat) const;

   /* Return point matrix of boundary element i of dimension Dim X #v, where for
      every vertex we give its coordinates in space of dimension Dim. */
   void GetBdrPointMatrix(int i, DenseMatrix &pointmat) const;

   /** This method aims to provide face information in a deciphered format, i.e.
       Mesh::FaceInformation, compared to the raw encoded information returned
       by Mesh::GetFaceElements() and Mesh::GetFaceInfos(). */
   FaceInformation GetFaceInformation(int f) const;

   void GetFaceElements (int Face, int *Elem1, int *Elem2) const;
   void GetFaceInfos (int Face, int *Inf1, int *Inf2) const;
   void GetFaceInfos (int Face, int *Inf1, int *Inf2, int *NCFace) const;

   /// @}

   /// @name Methods related to mesh partitioning
   /// @{

   /// @note The returned array should be deleted by the caller.
   int *CartesianPartitioning(int nxyz[]);
   /// @note The returned array should be deleted by the caller.
   int *GeneratePartitioning(int nparts, int part_method = 1);
   /// @todo This method needs a proper description
   void CheckPartitioning(int *partitioning_);

   /// @}

   /// Print the mesh to the given stream using the default MFEM mesh format.
   /// \see mfem::ofgzstream() for on-the-fly compression of ascii outputs. If
   /// @a comments is non-empty, it will be printed after the first line of the
   /// file, and each line should begin with '#'.
   virtual void Print(std::ostream &os = mfem::out,
                      const std::string &comments = "") const
   { Printer(os, "", comments); }

   /// Save the mesh to a file using Mesh::Print. The given @a precision will be
   /// used for ASCII output.
   virtual void Save(const std::string &fname, int precision=16) const;
};

/** Overload operator<< for std::ostream and Mesh; valid also for the derived
    class ParMesh */
std::ostream &operator<<(std::ostream &os, const Mesh &mesh);

/// @brief Print function for Mesh::FaceInformation.
std::ostream& operator<<(std::ostream &os, const Mesh::FaceInformation& info);


/** @brief Class containing a minimal description of a part (a subset of the
    elements) of a Mesh and its connectivity to other parts.

    The main purpose of this class is to facilitate the partitioning of serial
    meshes (in serial, i.e. on one processor) and save the parts in parallel
    MFEM mesh format.

    Another potential futrure purpose of this class could be to facilitate
    exchange of MeshParts between MPI ranks for repartitioning purposes. It can
    also potentially be used to implement parallel mesh I/O functions with
    partitionings that have number of parts different from the number of MPI
    tasks.

    @note Parts of NURBS or non-conforming meshes cannot be fully described by
    this class alone with its current data members. Such extensions may be added
    in the future.
*/
class MeshPart
{
protected:
   struct Entity { int geom; int num_verts; const int *verts; };
   struct EntityHelper
   {
      int dim, num_entities;
      int geom_offsets[Geometry::NumGeom+1];
      typedef const Array<int> entity_to_vertex_type[Geometry::NumGeom];
      entity_to_vertex_type &entity_to_vertex;

      EntityHelper(int dim_,
                   const Array<int> (&entity_to_vertex_)[Geometry::NumGeom]);
      Entity FindEntity(int bytype_entity_id);
   };

public:
   /// Reference space dimension of the elements
   int dimension;

   /// Dimension of the physical space into which the MeshPart is embedded.
   int space_dimension;

   /// Number of vertices
   int num_vertices;

   /// Number of elements with reference space dimension equal to 'dimension'.
   int num_elements;

   /** @brief Number of boundary elements with reference space dimension equal
       to 'dimension'-1. */
   int num_bdr_elements;

   /**
      Each 'entity_to_vertex[geom]' describes the entities of Geometry::Type
      'geom' in terms of their vertices. The number of entities of type 'geom'
      is:

          num_entities[geom] = size('entity_to_vertex[geom]')/num_vertices[geom]

      The number of all elements, 'num_elements', is:

          'num_elements' = sum_{dim[geom]=='dimension'} num_entities[geom]

      and the number of all boundary elements, 'num_bdr_elements' is:

          'num_bdr_elements' = sum_{dim[geom]=='dimension'-1} num_entities[geom]

      Note that 'entity_to_vertex' does NOT describe all "faces" in the mesh
      part (i.e. all 'dimension'-1 entities) but only the boundary elements.
      Also, note that lower dimesional entities ('dimension'-2 and lower) are
      NOT described by the respective array, i.e. the array will be empty.
   */
   Array<int> entity_to_vertex[Geometry::NumGeom];

   /** @brief Store the refinement flags for tetraheral elements. If all tets
       have zero refinement flags then this array is empty, i.e. has size 0. */
   Array<int> tet_refine_flags;

   /**
      Terminology: "by-type" element/boundary ordering: ordered by
      Geometry::Type and within each Geometry::Type 'geom' ordered as in
      'entity_to_vertex[geom]'.

      Optional re-ordering of the elements that will be used by (Par)Mesh
      objects constructed from this MeshPart. This array maps "natural" element
      ids (used by the Mesh/ParMesh objects) to "by-type" element ids (see
      above):

          "by-type" element id = element_map["natural" element id]

      The size of the array is either 'num_elements' or 0 when no re-ordering is
      needed (then "by-type" id == "natural" id).
   */
   Array<int> element_map;

   /// Optional re-ordering for the boundary elements, similar to 'element_map'.
   Array<int> boundary_map;

   /**
      Element attributes. Ordered using the "natural" element ordering defined
      by the array 'element_map'. The size of this array is 'num_elements'.
   */
   Array<int> attributes;

   /**
      Boundary element attributes. Ordered using the "natural" boundary element
      ordering defined by the array 'boundary_map'. The size of this array is
      'num_bdr_elements'.
   */
   Array<int> bdr_attributes;

   /**
      Optional vertex coordinates. The size of the array is either

          size = 'space_dimension' * 'num_vertices'

      or 0 when the vertex coordinates are not used, i.e. when the MeshPart uses
      a nodal GridFunction to describe its location in physical space. This
      array uses Ordering::byVDIM: "X0,Y0,Z0, X1,Y1,Z1, ...".
   */
   Array<real_t> vertex_coordinates;

   /**
      Optional serial Mesh object constructed on demand using the method
      GetMesh(). One use case for it is when one wants to construct FE spaces
      and GridFunction%s on the MeshPart for saving or MPI communication.
   */
   std::unique_ptr<Mesh> mesh;

   /** @name Connectivity to other MeshPart objects */
   ///@{

   /// Total number of MeshParts
   int num_parts;

   /** @brief Index of the part described by this MeshPart:
       0 <= 'my_part_id' < 'num_parts' */
   int my_part_id;

   /**
      A group G is a subset of the set { 0, 1, ..., 'num_parts'-1 } for which
      there is a mesh entity E (of any dimension) in the global mesh such that
      G is the set of the parts assigned (by the partitioning array) to the
      elements adjacent to E. The MeshPart describes only the "neighbor" groups,
      i.e. the groups that contain 'my_part_id'. The Table 'my_groups' defines
      the "neighbor" groups in terms of their part ids. In other words, it maps
      "neighbor" group ids to a (sorted) list of part ids. In particular, the
      number of "neighbor" groups is given by 'my_groups.Size()'. The "local"
      group { 'my_part_id' } has index 0 in 'my_groups'.
   */
   Table my_groups;

   /**
      Shared entities for this MeshPart are mesh entities of all dimensions less
      than 'dimension' that are generated by the elements of this MeshPart and
      at least one other MeshPart.

      The Table 'group_shared_entity_to_vertex[geom]' defines, for each group,
      the shared entities of Geometry::Type 'geom'. Each row (corresponding to a
      "neighbor" group, as defined by 'my_groups') in the Table defines the
      shared entities in a way similar to the arrays 'entity_to_vertex[geom]'.
      The "local" group (with index 0) does not have any shared entities, so the
      0-th row in the Table is always empty.

      IMPORTANT: the descriptions of the groups in this MeshPart must match
      their descriptions in all neighboring MeshParts. This includes the
      ordering of the shared entities within the group, as well as the vertex
      ordering of each shared entity.
   */
   Table group_shared_entity_to_vertex[Geometry::NumGeom];

   ///@}

   /** @brief Write the MeshPart to a stream using the parallel format
       "MFEM mesh v1.2". */
   void Print(std::ostream &os) const;

   /** @brief Construct a serial Mesh object from the MeshPart.

       The nodes of 'mesh' are NOT initialized by this method, however, the
       nodal FE space and nodal GridFunction can be created and then attached to
       the 'mesh'. The Mesh is constructed only if 'mesh' is empty, otherwise
       the method simply returns the object held by 'mesh'.
   */
   Mesh &GetMesh();
};


/** @brief Class that allows serial meshes to be partitioned into MeshPart
    objects, typically one MeshPart at a time, which can then be used to write
    the local mesh in parallel MFEM mesh format.

    Sample usage of this class: partition a serial mesh and save it in parallel
    MFEM format:
    \code
       // The array 'partitioning' can be obtained e.g. from
       // mesh->GeneratePartitioning():
       void usage1(Mesh *mesh, int num_parts, int *partitioning)
       {
          MeshPartitioner partitioner(*mesh, num_parts, partitioning);
          MeshPart mesh_part;
          for (int i = 0; i < num_parts; i++)
          {
             partitioner.ExtractPart(i, mesh_part);
             ofstream omesh(MakeParFilename("my-mesh.", i));
             mesh_part.Print(omesh);
          }
       }
    \endcode

    This class can also be used to partition a mesh and GridFunction(s) and save
    them in parallel:
    \code
       // The array 'partitioning' can be obtained e.g. from
       // mesh->GeneratePartitioning():
       void usage2(Mesh *mesh, int num_parts, int *partitioning,
                   GridFunction *gf)
       {
          MeshPartitioner partitioner(*mesh, num_parts, partitioning);
          MeshPart mesh_part;
          for (int i = 0; i < num_parts; i++)
          {
             partitioner.ExtractPart(i, mesh_part);
             ofstream omesh(MakeParFilename("my-mesh.", i));
             mesh_part.Print(omesh);
             auto lfes = partitioner.ExtractFESpace(mesh_part, *gf->FESpace());
             auto lgf = partitioner.ExtractGridFunction(mesh_part, *gf, *lfes);
             ofstream ofield(MakeParFilename("my-field.", i));
             lgf->Save(ofield);
          }
       }
    \endcode
*/
class MeshPartitioner
{
protected:
   Mesh &mesh;
   Array<int> partitioning;
   Table part_to_element;
   Table part_to_boundary;
   Table edge_to_element;
   Table vertex_to_element;

public:
   /** @brief Construct a MeshPartitioner.

       @param[in] mesh_         Mesh to be partitioned into MeshPart%s.
       @param[in] num_parts_    Number of parts to partition the mesh into.
       @param[in] partitioning_ Partitioning array: for every element in the
                                mesh gives the partition it belongs to; if NULL,
                                partitioning will be generated internally by
                                calling Mesh::GeneratePartitioning().
       @param[in] part_method   Partitioning method to be used in the call to
                                Mesh::GeneratePartitioning() when the provided
                                input partitioning is NULL.
   */
   MeshPartitioner(Mesh &mesh_, int num_parts_,
                   const int *partitioning_ = nullptr, int part_method = 1);

   /** @brief Construct a MeshPart corresponding to the given @a part_id.

       @param[in]  part_id    Partition index to extract; valid values are in
                              the range [0, num_parts).
       @param[out] mesh_part  Output MeshPart object; its contents is
                              overwritten, while potentially reusing existing
                              dynamic memory allocations.
   */
   void ExtractPart(int part_id, MeshPart &mesh_part) const;
};

// shift cyclically 3 integers left-to-right
inline void ShiftRight(int &a, int &b, int &c)
{
   int t = a;
   a = c;  c = b;  b = t;
}

}

#endif
