/*--------------------------------------------------------------------------------------+
|
|  Copyright: (c) 2011 Bentley Systems, Incorporated. All rights reserved. 
|
|  Redistribution and use in source and binary forms, with or without modification, 
|  are permitted provided that the following conditions are met:
|    
|    - Redistributions of source code must retain the above copyright notice, 
|       this list of conditions and the following disclaimer.
|    - Redistributions in binary form must reproduce the above copyright notice, 
|       this list of conditions and the following disclaimer in the documentation 
|       and/or other materials provided with the distribution.
|    - Neither the name of Bentley Systems, Incorporated nor the names of its 
|       contributors may be used to endorse or promote products derived from this 
|       software without specific prior written permission.
|
|  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
|  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
|  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
|  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
|  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
|  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
|  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
|  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN 
|  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH 
|  DAMAGE.
|
+--------------------------------------------------------------------------------------*/

using System;
using System.Diagnostics;
using System.Collections.Generic;
using Bentley.Geometry;
using Bentley.GenerativeComponents;
using Bentley.GenerativeComponents.GeneralPurpose;
using Bentley.GenerativeComponents.GCScript;
using Bentley.GenerativeComponents.MicroStation;
using Bentley.Interop.MicroStationDGN;
using Bentley.GenerativeComponents.Features.Specific;

namespace Bentley.GenerativeComponents.Features.Specific
{
    [GCNamespace("User")]
    public class MyFabricationPlanning : PolygonalGeometry, IPolygon
    {
        private double m_lenght = 0;
        private double m_area = 0;

        public MyFabricationPlanning(): base(){}
        public MyFabricationPlanning(Feature parentFeature) : base(parentFeature) { }

        protected  override Feature NewInstance() { return new MyFabricationPlanning(this); }

        //public override DPoint3d[] IPolygon_Vertices
        //{
        //    get
        //    {

        //        DPoint3d[] returnVertices = null;

        //        if ((InternalVertices != null) && (InternalVertices.Length > 0))
        //        {
        //            int iCount = InternalVertices.Length;

        //            returnVertices = new DPoint3d[iCount];

        //            for (int i = 0; i < iCount; i++)
        //            {
        //                returnVertices[i] = new DPoint3d(InternalVertices[i].X, InternalVertices[i].Y, InternalVertices[i].Z);
        //            }
        //        }

        //        return returnVertices;
        //    }
        //}
        DPoint3d[] IPolygon.Vertices { get { return IPolygon_Vertices; } }


        ///<summary>Layout an array of Curves in the XYplane of a specified CoordinateSystem</summary>
        ///<categories>Curve Based</categories>
        ///<param name="CoordinateSystem">CoordinateSystem defining layout space</param>
        ///<param name="CurveArray">Array of Curves to layout</param>
        ///<param name="MaxRowLength">Maximum number of Curves in one Layout Row</param>
        ///<param name="Xspacing">Spacing in Xdirection between Curve Centroids</param>
        ///<param name="Yspacing">Spacing in Ydirection between Curve Centroids</param>
        ///<param name="PlaneArray">Optional List of Planes corresponding to Curve List ensures layout Curves are co-planar with XY plane of CoordinateSystem</param>
        [Update]
        public virtual bool LayoutCurveArray
        (
        FeatureUpdateContext updateContext,
        [ParentModel]  CoordinateSystem CoordinateSystem,
        IGCObject  CurveArray,
        // [Optional]  int      CurveToUseToDefinedRotation,
        // [Optional]double   RotationAngle,
        int     MaxRowLength,
        double  Xspacing,
        double  Yspacing,
        [Optional]IPlane[]  PlaneArray
        )
        {
            if((CoordinateSystem!=null)&&(CurveArray!=null)&&(CurveArray.IsList((Type) null))&&(CurveArray.Count((Type) null)>0))
            {

                Point3d com_cellOrigin = DgnTools.ToPoint3d(DPoint3d.Zero);

                DMatrix3d layoutOrientation = CoordinateSystem.DTransform3d.Matrix;

                DPoint3d origin = CoordinateSystem.DPoint3d;

                int iCount = CurveArray.Count((Type) null);

                int iPlaneCount = 0;

                bool bUsePlanes = false;

                if(PlaneArray!=null)
                {
                    iPlaneCount = PlaneArray.Length;

                    if(iPlaneCount==iCount)
                    {
                        bUsePlanes = true;
                    }
                }


                /*
                int iCurveToUse = CurveToUseToDefinedRotation;

                if(CurveToUseToDefinedRotation<0) iCurveToUse = 0;

                if(CurveToUseToDefinedRotation>(iCount-1)) iCurveToUse = iCount;
                */

                int iMaxRowLengthToUse = MaxRowLength;

                if(MaxRowLength<1) iMaxRowLengthToUse = 1;

                int iRows = iCount/iMaxRowLengthToUse;

                int iRemainder = iCount-(MaxRowLength*iRows);

                if(iRemainder>0) iRows++;

                Element[] rowCell = new Element[iRows];

                // correct origin

                origin.Y = origin.Y + iRows * Yspacing;

                double originalX = origin.X;

                /*
                Transform3d inputTransform = CurveArray[iCurveToUse].m_DTransform3d;

                Transform3d inverseTransform = GeometryTools.m_DTransform3dInverse( inputTransform);
                */

                DTransform3d layoutTransform = DTransform3d.Identity;

                layoutTransform.InitMatrixAndTranslation(ref layoutOrientation, ref origin);

                Transform3d com_layoutTransform = DgnTools.ToTransform3d(layoutTransform);

                int iCurveIndex = 0;

                m_lenght = m_area = 0;//TR276745
                for(int i=0; i<iRows; i++)
                {
                    Element[] rowElement = new Element[iMaxRowLengthToUse];

                    for(int j=0; j<iMaxRowLengthToUse; j++)
                    {
                        ICurve nestedFeature = CurveArray.GetItem((Type) null, iCurveIndex) as ICurve;

                        if((nestedFeature.IsList((Type) null))&&(nestedFeature.Count((Type) null)>1))
                        {
                            int kCount = nestedFeature.Count((Type) null);

                            Element[] nestedNestedElement = new Element[kCount];

                            for (int k=0; k<kCount; k++)
                            {
                                ICurve nestedNestedFeature = nestedFeature.GetItem((Type) null, k) as ICurve;

                                //TR276745
                                m_lenght += nestedNestedFeature.Length;
                                m_area += nestedNestedFeature.Area;

                                if(nestedNestedFeature!=null)
                                {
                                    nestedNestedElement[k] = DgnTools.CloneElement(nestedNestedFeature.GetElement());

                                    DTransform3d inputTransform = ((Feature)nestedNestedFeature).m_DTransform3d;

                                    if(bUsePlanes)
                                    {
                                        inputTransform = ((Feature)PlaneArray[iCurveIndex]).m_DTransform3d;
                                    }

                                    DTransform3d inverseTransform = GeometryTools.InverseOf(inputTransform);

                                    Transform3d com_inverseTransform = DgnTools.ToTransform3d(inverseTransform);

                                    nestedNestedElement[k].IsLocked = false;
                                    nestedNestedElement[k].Transform(ref com_inverseTransform);
                                }
                            }

                            rowElement[j] = (Element) DgnTools.CreateCellElement1("nested[" + j +"]", nestedNestedElement, COM_ZERO_POINT3D, false);
                        }
                        else
                        {
                            rowElement[j] =  DgnTools.CloneElement(nestedFeature.GetElement());

                            DTransform3d inputTransform = ((Feature)nestedFeature).m_DTransform3d;

                            if(bUsePlanes)
                            {
                                inputTransform = ((Feature)PlaneArray[iCurveIndex]).m_DTransform3d;
                            }

                            DTransform3d inverseTransform = GeometryTools.InverseOf(inputTransform);

                            Transform3d com_inverseTransform = DgnTools.ToTransform3d(inverseTransform);

                            rowElement[j].IsLocked = false;
                            rowElement[j].Transform(ref com_inverseTransform);
                        }

                        com_layoutTransform = DgnTools.ToTransform3d(layoutTransform);

                        rowElement[j].Transform(ref com_layoutTransform);

                        if(iCurveIndex==iCount-1) break;

                        iCurveIndex++;

                        origin.X = origin.X + Xspacing;

                        layoutTransform.InitMatrixAndTranslation(ref layoutOrientation, ref origin);
                    }

                    origin.Y = origin.Y - Yspacing;

                    origin.X = originalX;

                    layoutTransform.InitMatrixAndTranslation(ref layoutOrientation, ref origin);

                    rowCell[i] = (Element) DgnTools.CreateCellElement1("row[" + i +"]", rowElement, COM_ZERO_POINT3D, false);
                }

                SetElement(DgnTools.CreateCellElement1("curve_array", rowCell, COM_ZERO_POINT3D, false));

                return true;
            }

            return false;
        }


        ///<summary>Unfold a rectangular array of Polygons into strips made planar by triangulation of vertices</summary>
        ///<categories>Polygon Based</categories>
        ///<param name="CoordinateSystem">CoordinateSystem defining layout space</param>
        ///<param name="Polygons">Rectangular Array of Polygons to layout</param>
        ///<param name="BoundaryEdgeColor">Color of Outer Edges (cutting line on laser cutter)</param>
        ///<param name="InternalEdgeColor">Color of Inner Edges (scoring line on laser cutter)</param>
        ///<param name="InterRowDistance">Distance between each Strip</param>
        /// <param name="CorrectForNonPlanarQuad">Unfold NonPlanar Polygons as Planar or NonPlanar</param>
        ///<param name="TextStyle">Assign optional TextStyle to format Polygon ID tags when Construction is Displayed</param>
        [Update]
        public virtual bool UnfoldPolygonsIntoPlanarStrips // previously LayoutPolygonsUsingVertexMapping
        (
        FeatureUpdateContext    updateContext,
        [ParentModel]           CoordinateSystem CoordinateSystem,
        [ArrayBuilder(2,99)]    Polygon[][] Polygons,
        int                     BoundaryEdgeColor,
        int                     InternalEdgeColor,
        double                  InterRowDistance,
        //[Optional]int[] IndicesOfVerticesDefiningPolygonPosition,
        //[Optional]int[] IndicesOfVerticesDefiningContextToPlaceNextPolygon,
        [Optional, DefaultValue(true)]bool  CorrectForNonPlanarQuad,
        [Optional] Bentley.GenerativeComponents.Features.Specific.TextStyle TextStyle
        )
        {
            if((CoordinateSystem!=null)&&(Polygons!=null)&&(Polygons.Length>0)&&(Polygons[0]!=null)&&(Polygons[0].Length>0))
            {
                //added for TR282274
                DeleteAllSubFeatures(updateContext);
                DeleteElement(updateContext);

                //added for TR278591
                //bool isQuads = (Polygons[0][0].Facet == FacetOption.Quads ? true : false);
                bool isQuads = true;//ww TBD : the Facet is internal property.

                int[] dimension = new int[Polygons.Length];
                for (int i = 0; i < Polygons.Length; i++)
                {
                    if (!CheckPolygonData(Polygons[i], ref dimension[i])) return false;
                    if (dimension[i] != 3 && dimension[i] != 4) return false;
                }

                //if (IndicesOfVerticesDefiningPolygonPosition == null) IndicesOfVerticesDefiningPolygonPosition = new int[3] { 0, 1, 3 };
                //if (IndicesOfVerticesDefiningContextToPlaceNextPolygon == null) IndicesOfVerticesDefiningContextToPlaceNextPolygon = new int[3] { 3, 2, 1 };
                int[] IndicesOfVerticesDefiningPolygonPosition = new int[3] { 0, 1, 3 };
                int[] IndicesOfVerticesDefiningContextToPlaceNextPolygon = new int[3] { 3, 2, 1 };

                SymbologyAndLevelUsage = SymbologyAndLevelUsageOption.AssignToElement;

                DMatrix3d layoutOrientation = CoordinateSystem.DTransform3d.Matrix;

                m_DTransform3d = CoordinateSystem.m_DTransform3d;

                DPoint3d origin = DPoint3d.Zero;

                origin.ProductOf(ref m_DTransform3d, ref origin);

                int iCount = Polygons.Length;

                Element[] lineElement = new Element[iCount];

                double rowPosition = 0.0;

                for(int i=0; i<iCount; i++)
                {
                    if (isQuads)
                    {
                        lineElement[i] = LayoutPolygonArrayIntoPlanarStrips(layoutOrientation,
                                         origin,
                                         Polygons[i],
                                         IndicesOfVerticesDefiningPolygonPosition,
                                         IndicesOfVerticesDefiningContextToPlaceNextPolygon,
                                         BoundaryEdgeColor,
                                         InternalEdgeColor,
                                         CorrectForNonPlanarQuad,
                                         true,
                                         TextStyle,
                                         "[" + i + "]");
                    }
                    else
                    {
                        //TR278591
                        lineElement[i] = LayoutPolygonsIntoPlanarStrips(layoutOrientation,
                                                          origin,
                                                          Polygons[i],
                                                          dimension[i],
                                                          BoundaryEdgeColor,
                                                          InternalEdgeColor,
                                                          CorrectForNonPlanarQuad,
                                                          true,
                                                          TextStyle,
                                                          "[" + i + "]");
                    }

                    rowPosition = rowPosition + InterRowDistance;

                    origin.Init(rowPosition, 0.0, 0.0);

                    origin.ProductOf(ref m_DTransform3d, ref origin);
                }

                SetElement(DgnTools.CreateCellElement1(Name, lineElement, COM_ZERO_POINT3D, false));
                return true;
            }
            return false;
        }
        private bool CheckPolygonData(Polygon[] polygons, ref int dimPolygon)
        {
            bool ok = true;

            dimPolygon = polygons[0].Vertices.Length;
            int iCount = polygons.Length;
            for (int i = 0; i < iCount; i++)
            {
                int number = polygons[i].Vertices.Length;
                if (dimPolygon != number)
                {
                    ok = false;
                    break;
                }
            }

            return ok;
        }

        private DPoint3d[] GetPolygonPoints(Polygon p, bool closed)
        {
            //int initialCount = PolygonArray[i].IPolygon_Vertices.Length; 
            int initialCount = p.Vertices.Length; //ww : TBD : 

            if (closed)
            {
                DPoint3d[] vertices = new DPoint3d[initialCount + 1];

                DPoint3d[] verticesToGet = p.VerticesAsDPoint3ds; // ww TBD : 
                for (int j = 0; j < initialCount; ++j)
                {
                    vertices[j] = verticesToGet[j];
                }

                vertices[initialCount] = verticesToGet[0];

                return vertices;
            }
            else
            {
                DPoint3d[] vertices = new DPoint3d[initialCount];

                DPoint3d[] verticesToGet = p.VerticesAsDPoint3ds; // ww TBD : 
                for (int j = 0; j < initialCount; ++j)
                {
                    vertices[j] = verticesToGet[j];
                }

                return vertices;
            }
        }


        internal Element LayoutPolygonArrayIntoPlanarStrips
        (
        DMatrix3d layoutOrientation,
        DPoint3d origin,
        Polygon[] PolygonArray,
        int[] inputIndices,
        int[] outputIndices,
        int boundaryEdgeColor,
        int internalEdgeColor,
        bool bCorrectForNonPlanarQuad,
        bool bHanding,
        Bentley.GenerativeComponents.Features.Specific.TextStyle TextStyle,
        string Key
        )
        {
            if ((PolygonArray != null) && (PolygonArray.Length > 0) && (inputIndices.Length == 3) && (outputIndices.Length == 3))
            {
                DVector3d ZVector = new DVector3d(0.0, 0.0, 1.0);

                int iCount = PolygonArray.Length;

                DTransform3d layoutTransform = DTransform3d.Identity;
                layoutTransform.InitMatrixAndTranslation(ref layoutOrientation, ref origin);

                int iResultCount = 0;

                int noOfLines = iCount * 4 + 1;

                if (bCorrectForNonPlanarQuad) noOfLines = noOfLines + iCount;

                Point3d com_start = DgnTools.ToPoint3d(DPoint3d.Zero);
                Point3d com_end = DgnTools.ToPoint3d(DPoint3d.Zero);

                Element[] lineElement = new Element[noOfLines];

                for (int i = 0; i < iCount; ++i)
                {
                    DPoint3d[] vertices = GetPolygonPoints(PolygonArray[i], true);

                    if ((vertices != null) && (vertices.Length > 3))
                    {
                        int iVertexCount = vertices.Length;

                        DMatrix3d inputOrientation = DMatrix3d.Identity;

                        inputOrientation.InitOrthogonalFrameFromOriginXY(ref vertices[inputIndices[0]], ref vertices[inputIndices[1]], ref vertices[inputIndices[2]]);

                        DTransform3d inputTransform = DTransform3d.Identity;

                        inputTransform.InitMatrixAndTranslation(ref inputOrientation, ref vertices[inputIndices[0]]);

                        DTransform3d inverseTransform = GeometryTools.InverseOf(inputTransform);

                        for (int j = 0; j < iVertexCount; j++)
                        {
                            vertices[j].ProductOf(ref inverseTransform, ref vertices[j]);

                            if (outputIndices[0] == j)
                            {
                                vertices[j].ProductOf(ref layoutTransform, ref vertices[j]);
                            }
                            else
                            {
                                vertices[j].ProductOf(ref layoutTransform, ref vertices[j]);
                            }
                        }

                        // -------------------

                        if (ConstructionsVisible == true)
                        {
                            if (!Key.IsEmpty())
                            {
                                Text Text = new Text(this);

                                Text.FromDTransform3d(layoutTransform, Key + "[" + i + "]", TextStyle);

                                AddConstructionFeature(Text, "Key");
                            }
                        }

                        // -------------------

                        if (bCorrectForNonPlanarQuad)
                        {
                            DMatrix3d correctionInputOrientation = DMatrix3d.Identity;
                            correctionInputOrientation.InitOrthogonalFrameFromOriginXY(ref vertices[outputIndices[0]], ref vertices[inputIndices[1]], ref vertices[outputIndices[1]]);

                            DTransform3d correctionInputTransform = DTransform3d.Identity;
                            correctionInputTransform.InitMatrixAndTranslation(ref correctionInputOrientation, ref vertices[outputIndices[0]]);

                            DTransform3d inverseCorrectionTransform = GeometryTools.InverseOf(correctionInputTransform);

                            DMatrix3d correctionOrientation = GeometryTools.GetMatrixFromPointsAndZVector(vertices[outputIndices[0]], vertices[inputIndices[1]], ZVector);

                            DTransform3d correctionTransform = DTransform3d.Identity;
                            correctionTransform.InitMatrixAndTranslation(ref correctionOrientation, ref vertices[outputIndices[0]]);

                            vertices[outputIndices[1]].ProductOf(ref inverseCorrectionTransform, ref vertices[outputIndices[1]]);
                            vertices[outputIndices[1]].ProductOf(ref correctionTransform, ref vertices[outputIndices[1]]);
                        }

                        if ((inputIndices != null) && (inputIndices.Length == 3))
                        {
                            for (int j = 1; j < iVertexCount; j++)
                            {
                                com_start = DgnTools.ToPoint3d(vertices[j - 1]);
                                com_end = DgnTools.ToPoint3d(vertices[j]);

                                lineElement[iResultCount] = (Element)MSApp.CreateLineElement2(TemplateElement, ref com_start, ref com_end);


                                if (((j == inputIndices[0]) && ((j - 1) == inputIndices[1]))
                                   || (((j - 1) == inputIndices[0]) && (j == inputIndices[1])))
                                {
                                    if (i == 0)                                                  // first tiime draw the first edge
                                    {
                                        lineElement[iResultCount].Color = boundaryEdgeColor;
                                        iResultCount++;
                                    }
                                    else
                                    {
                                        lineElement[iResultCount].Color = internalEdgeColor;
                                        iResultCount++;
                                    }

                                }

                                else if (((j == outputIndices[0]) && ((j - 1) == outputIndices[1]))
                                        || (((j - 1) == outputIndices[0]) && (j == outputIndices[1])))
                                {
                                    if (i == iCount - 1)                                           // last tiime draw the last edge
                                    {
                                        lineElement[iResultCount].Color = boundaryEdgeColor;
                                        iResultCount++;
                                    }
                                    else
                                    {
                                        // do nothing
                                    }
                                }

                                else
                                {
                                    lineElement[iResultCount].Color = boundaryEdgeColor;
                                    iResultCount++;
                                }
                            }

                            if (bCorrectForNonPlanarQuad)
                            {
                                if (bHanding)
                                {
                                    com_start = DgnTools.ToPoint3d(vertices[inputIndices[1]]);
                                    com_end = DgnTools.ToPoint3d(vertices[outputIndices[0]]);

                                    lineElement[iResultCount] = (Element)MSApp.CreateLineElement2(TemplateElement, ref com_start, ref com_end);
                                }

                                else
                                {
                                    com_start = DgnTools.ToPoint3d(vertices[inputIndices[0]]);
                                    com_end = DgnTools.ToPoint3d(vertices[outputIndices[1]]);

                                    lineElement[iResultCount] = (Element)MSApp.CreateLineElement2(TemplateElement, ref com_start, ref com_end);
                                }

                                lineElement[iResultCount].Color = internalEdgeColor;
                                iResultCount++;
                            }
                        }
                    }

                    DMatrix3d outputOrientation = GeometryTools.GetMatrixFromPointsAndZVector(vertices[outputIndices[0]], vertices[outputIndices[1]], ZVector);

                    layoutTransform.InitMatrixAndTranslation(ref outputOrientation, ref vertices[outputIndices[0]]);
                }

                com_start = DgnTools.ToPoint3d(origin);

                return (Element)DgnTools.CreateCellElement1(Name, lineElement, com_start, false);
            }
            return null;
        }
        internal Element LayoutPolygonsIntoPlanarStrips
        (
        DMatrix3d layoutOrientation,
        DPoint3d origin,
        Polygon[] PolygonArray,
        int dimension,
        int boundaryEdgeColor,
        int internalEdgeColor,
        bool bCorrectForNonPlanarQuad,
        bool bHanding,
        Bentley.GenerativeComponents.Features.Specific.TextStyle TextStyle,
        string Key
        )
        {
            DVector3d ZVector = new DVector3d(0.0, 0.0, 1.0);

            DTransform3d layoutTransform = DTransform3d.Identity;
            layoutTransform.InitMatrixAndTranslation(ref layoutOrientation, ref origin);

            int numtotalpoints = PolygonArray.Length * dimension;
            DPoint3d[] totalpoints = new DPoint3d[numtotalpoints];
            for (int i = 0; i < PolygonArray.Length; ++i)
            {
                DPoint3d[] vertices = PolygonArray[i].VerticesAsDPoint3ds;
                for (int j = 0; j < dimension; ++j)
                {
                    totalpoints[i * dimension + j] = vertices[j];
                }
            }

            DTransform3d inputTransform = DTransform3d.Identity;
            DVector3d vector = new DVector3d(ref origin, ref totalpoints[0]);
            inputTransform.InitTranslation(ref vector);
            DTransform3d inverseTransform = GeometryTools.InverseOf(inputTransform);

            for (int i = 0; i < numtotalpoints; i++)
            {
                totalpoints[i].ProductOf(ref inverseTransform, ref totalpoints[i]);
            }

            Point3d com_start = DgnTools.ToPoint3d(DPoint3d.Zero);
            Point3d com_end = DgnTools.ToPoint3d(DPoint3d.Zero);

            int noOfLines = PolygonArray.Length * dimension;
            if (dimension == 4) noOfLines = noOfLines + PolygonArray.Length;
            //if (bCorrectForNonPlanarQuad) noOfLines = noOfLines + PolygonArray.Length;
            Element[] lineElement = new Element[noOfLines];

            int iResultCount = 0;
            for (int i = 0; i < PolygonArray.Length; i++)
            {

                DPoint3d[] vertices = new DPoint3d[dimension + 1];
                for (int j = 0; j < dimension; j++)
                {
                    vertices[j] = totalpoints[i * dimension + j];
                }
                vertices[dimension] = vertices[0];

                if (ConstructionsVisible == true)
                {
                    if (!Key.IsEmpty())
                    {
                        Text Text = new Text(this);

                        Text.FromDTransform3d(layoutTransform, Key + "[" + i + "]", TextStyle);

                        AddConstructionFeature(Text, "Key");
                    }
                }


                if (bCorrectForNonPlanarQuad)
                {
                    for (int j = 0; j < dimension + 1; j++)
                    {
                        DVector3d vectorz = new DVector3d(0.0, 0.0, vertices[j].Z);
                        DTransform3d inputtransform = DTransform3d.Identity;
                        inputTransform.InitTranslation(ref vectorz);
                        vertices[j].ProductOfOrthogonalInverse(ref inputTransform, ref vertices[j]);
                    }
                }

                if (vertices == null || vertices.Length < 3) return null;

                for (int j = 1; j < dimension + 1; j++)
                {
                    com_start = DgnTools.ToPoint3d(vertices[j - 1]);
                    com_end = DgnTools.ToPoint3d(vertices[j]);

                    lineElement[iResultCount] = (Element)MSApp.CreateLineElement2(TemplateElement, ref com_start, ref com_end);
                    lineElement[iResultCount].Color = boundaryEdgeColor;
                    iResultCount++;
                }

                if (dimension == 4)
                {
                    com_start = DgnTools.ToPoint3d(vertices[0]);
                    com_end = DgnTools.ToPoint3d(vertices[2]);
                    lineElement[iResultCount] = (Element)MSApp.CreateLineElement2(TemplateElement, ref com_start, ref com_end);
                    lineElement[iResultCount].Color = internalEdgeColor;
                    iResultCount++;
                }

                DMatrix3d outputOrientation = GeometryTools.GetMatrixFromPointsAndZVector(vertices[0], vertices[1], ZVector);

                layoutTransform.InitMatrixAndTranslation(ref outputOrientation, ref vertices[0]);
            }

            com_start = DgnTools.ToPoint3d(origin);

            return (Element)DgnTools.CreateCellElement1(Name, lineElement, com_start, false);
        }

        ///<summary>Layout an array of Polygons in the XYplane of a specified CoordinateSystem</summary>
        ///<categories>Polygon Based</categories>
        ///<param name="CoordinateSystem">CoordinateSystem defining Layout Space</param>
        ///<param name="Polygons">Rectangular configuration of Polygons to Layout</param>
        ///<param name="Xspacing">Spacing in Xdirection between Polygon Centroids</param>
        ///<param name="Yspacing">Spacing in Ydirection between Polygon Centroids</param>
        ///<param name="TextStyle">Assign optional TextStyle to format Polygon ID tags when Construction is Displayed</param>
        [Update]
        public virtual bool LayoutPolygons
        (
        FeatureUpdateContext updateContext,
        [ParentModel]  CoordinateSystem CoordinateSystem,
        [ArrayBuilder(2,99)]Polygon[][] Polygons,
        double  Xspacing,
        double  Yspacing,
        [Optional] bool Fill,
        // [Optional] bool ForcePlanar,
        [Optional] Bentley.GenerativeComponents.Features.Specific.TextStyle TextStyle,
        [Optional, DefaultValue(false)] bool DisableOrientaion
        )
        {
            if((CoordinateSystem!=null)&&(Polygons!=null)&&(Polygons.Length>0)&&(Polygons[0]!=null)&&(Polygons[0].Length>0))
            {
                //added for TR282274
                DeleteAllSubFeatures(updateContext);
                DeleteElement(updateContext);

                DeleteConstituentFeatures(updateContext);

                DMatrix3d layoutOrientation = CoordinateSystem.DTransform3d.Matrix;

                DPoint3d origin = CoordinateSystem.DPoint3d;

                DTransform3d layoutTransform = DTransform3d.Identity;

                int iCount = Polygons.Length;

                double Xinc = 0.0; double Yinc = 0.0;

                double Xmin=0.0; double Ymin=0.0; double Xmax=0.0; double Ymax=0.0;

                FindMaxRangeOfPolygonInPolygonArray( ref Xmin, ref Xmax,ref Ymin,ref Ymax, Polygons[0], true);

                for(int i=1; i<iCount; i++)
                {
                    FindMaxRangeOfPolygonInPolygonArray( ref Xmin, ref Xmax,ref Ymin,ref Ymax, Polygons[i], false);
                }

                Xinc = Xmax - Xmin;
                Yinc = Ymax - Ymin;

                if(Xspacing>0.0001) Xinc = Xinc + Xspacing;
                if(Yspacing>0.0001) Yinc = Yinc + Yspacing;

                LetConstituentFeaturesBeDirectlyIndexible();

                for(int i=0; i<iCount; i++)
                {
                    MyFabricationPlanning nestedFabricationPlanning = new MyFabricationPlanning(this);

                    nestedFabricationPlanning.AlignOptions(this);
                    if (nestedFabricationPlanning.LayoutPolygonArray(updateContext, Xinc, Yinc, layoutOrientation, origin, Polygons[i], TextStyle, "[" + i + "]", DisableOrientaion))
                    {
                        nestedFabricationPlanning.DenoteAsHavingBeenUpdated();
                        origin.Init(origin.X + Xinc, origin.Y, origin.Z);

                        nestedFabricationPlanning.SetSuccess(true);
                        m_area += calcAreaOfPolygons(Polygons[i]);
                        SetConstituentFeature(updateContext, i, nestedFabricationPlanning);
                    }
                }
                return true;
            }
            return false;
        }

        private double calcAreaOfPolygons(Polygon[] polygons)
        {
            double total = 0;
            int nCount = polygons.Length;
            for (int i=0; i<nCount; i++)
            {
                total += polygons[i].Area;
            }
            return total;
        }

        private bool calcAreaAndLength(out double length, out double area)
        {
            length = 0;
            area = 0;
            Element elem = GetElement();
            if (elem != null)
            {
                if (elem.IsShapeElement())
                {
                    ShapeElement shape = elem.AsShapeElement();
                    area = shape.Area();
                    length = shape.Perimeter();
                    return true;
                }
                else
                {
                    if (elem.IsCurveElement())
                    {
                        CurveElement curve =elem.AsCurveElement();
                        area = 0;
                        length = curve.Length;
                        return true;
                    }
                }
            }
            return false;
        }

        [PropertyCalculator("Length")]
        public override void LengthPropertyCalculator ()
        {
            calcAreaAndLength(out m_lenght, out m_area);
            Length = m_lenght;
        }
        [TotalPropertyCalculator("TotalLength")]
        public override void TotalLengthPropertyCalculator ()
        {
            TotalLength = Length;
        }
        [PropertyCalculator("Area")]
        public override void AreaPropertyCalculator ()
        {
            calcAreaAndLength(out m_lenght, out m_area);
            Area = m_area;
        }
        [TotalPropertyCalculator("TotalArea")]
        public override void TotalAreaPropertyCalculator ()
        {
            TotalArea = Area;
        }

        internal void FindMaxRangeOfPolygonInPolygonArray
        (
        ref double Xmin,
        ref double Xmax,
        ref double Ymin,
        ref double Ymax,
        Polygon[] PolygonArray,
        bool bFirstTime
        )
        {

            Xmin = 0.0; Ymin = 0.0; Xmax = 0.0; Ymax = 0.0;

            if ((PolygonArray != null) && (PolygonArray.Length > 0))
            {
                int iCount = PolygonArray.Length;

                for (int i = 0; i < iCount; i++)
                {
                    DPoint3d[] vertices = GetPolygonPoints (PolygonArray[i], false);

                    if ((vertices != null) && (vertices.Length > 2))
                    {
                        int iVertexCount = vertices.Length;

                        DMatrix3d orientation = DMatrix3d.Identity;

                        orientation.InitOrthogonalFrameFromOriginXY(ref vertices[0], ref vertices[1], ref vertices[2]);

                        DTransform3d localTransform = DTransform3d.Identity;

                        localTransform.InitMatrixAndTranslation(ref orientation, ref vertices[0]);

                        DTransform3d inverseTransform = GeometryTools.InverseOf(localTransform);

                        for (int j = 0; j < iVertexCount; j++)
                        {
                            vertices[j].ProductOf(ref inverseTransform, ref vertices[j]);

                            if ((i == 0) && (j == 0))
                            {
                                Xmin = vertices[j].X;
                                Ymin = vertices[j].Y;
                                Xmax = vertices[j].X;
                                Ymax = vertices[j].Y;
                            }
                            else
                            {
                                if (Xmin > vertices[j].X) Xmin = vertices[j].X;
                                if (Ymin > vertices[j].Y) Ymin = vertices[j].Y;
                                if (Xmax < vertices[j].X) Xmax = vertices[j].X;
                                if (Ymax < vertices[j].Y) Ymax = vertices[j].Y;
                            }
                        }
                    }
                }
            }
        }
        internal bool LayoutPolygonArray
        (
        FeatureUpdateContext updateContext,
        double Xinc,
        double Yinc,
        DMatrix3d layoutOrientation,
        DPoint3d origin,
        Polygon[] PolygonArray,
        Bentley.GenerativeComponents.Features.Specific.TextStyle TextStyle,
        string Key,
        bool disableOrientation
        )
        {
            if ((PolygonArray != null) && (PolygonArray.Length > 0))
            {
                DeleteConstituentFeatures(updateContext);
                LetConstituentFeaturesBeDirectlyIndexible();

                DTransform3d layoutTransform = DTransform3d.Identity;

                int iCount = PolygonArray.Length;
                for (int i = 0; i < iCount; i++)
                {
                    DPoint3d origin_tmp = new DPoint3d(origin.X, origin.Y + Yinc * i, origin.Z);
                    layoutTransform.InitMatrixAndTranslation(ref layoutOrientation, ref origin_tmp);

                    DPoint3d[] vertices = GetPolygonPoints(PolygonArray[i], false);
                    int iVertexCount = vertices.Length;

                    DTransform3d localTransform = DTransform3d.Identity;
                    DMatrix3d orientation = DMatrix3d.Identity;
                    if (!disableOrientation)
                    {
                        if (iVertexCount > 3)
                            orientation.InitOrthogonalFrameFromOriginXY(ref vertices[0], ref vertices[1], ref vertices[3]);
                        else if (iVertexCount > 2)
                            orientation.InitOrthogonalFrameFromOriginXY(ref vertices[0], ref vertices[1], ref vertices[2]);
                    }
                    localTransform.InitMatrixAndTranslation(ref orientation, ref vertices[0]);
                    DTransform3d inverseTransform = GeometryTools.InverseOf(localTransform);

                    for (int j = 0; j < iVertexCount; j++)
                    {
                        vertices[j].ProductOf(ref inverseTransform, ref vertices[j]);
                        vertices[j].ProductOf(ref layoutTransform, ref vertices[j]);
                    }

                    MyFabricationPlanning nestedFabricationPlanning = new MyFabricationPlanning(this);

                    if (nestedFabricationPlanning.ComputeByVertices(GetStandardFeatureUpdateContext(), vertices, false, TextStyle, Key + "[" + i + "]"))
                    {
                        nestedFabricationPlanning.AlignOptions(this);
                        nestedFabricationPlanning.DenoteAsHavingBeenUpdated();
                        nestedFabricationPlanning.SetSuccess(true);
                        SetConstituentFeature(updateContext, i, nestedFabricationPlanning);
                    }
                }

                return true;
            }
            return false;
        }



        ///<summary>Layout an array of Polygons in the XYplane of a specified CoordinateSystem</summary>
        ///<categories>Polygon Based</categories>
        ///<param name="CoordinateSystem">CoordinateSystem defining layout space</param>
        ///<param name="PolygonArray">Rectangular array of Polygons to Layout</param>
        ///<param name="Xspacing">Spacing in Xdirection between Polygon Centroids</param>
        ///<param name="Yspacing">Spacing in Ydirection between Polygon Centroids</param>
        ///<param name="TextStyle">Assign optional TextStyle to format Polygon ID tags when Construction is Displayed</param>
        // [Update]
        public virtual bool LayoutPolygonArray
        (
        FeatureUpdateContext updateContext,
        [ParentModel]  CoordinateSystem CoordinateSystem,
        [ArrayBuilder(2,99)]Polygon[] PolygonArray,
        double  Xspacing,
        double  Yspacing,
        [Optional] bool Fill,
        [Optional]TextStyle TextStyle
        )
        {
            //if((CoordinateSystem!=null)&&(PolygonArray!=null)&&(PolygonArray.Length>0))
            //{
            //    double Xmin=0.0; double Ymin=0.0; double Xmax=0.0; double Ymax=0.0;

            //    FindMaxRangeOfPolygonInPolygonArray( ref Xmin, ref Xmax,ref Ymin,ref Ymax, PolygonArray,true);

            //    double Xinc = 0.0;
            //    double Yinc = 0.0;

            //    if(Xspacing>0.0001) Xinc = Xmax - Xmin + Xspacing;
            //    if(Yspacing>0.0001) Yinc = Ymax - Ymin + Yspacing;

            //    DMatrix3d layoutOrientation = CoordinateSystem.DTransform3d.Matrix;

            //    DPoint3d origin = CoordinateSystem.DPoint3d;

            //    return LayoutPolygonArray(updateContext, Xinc, Yinc, layoutOrientation, origin, PolygonArray, TextStyle, "", DisableOrientaion);
            //}
            return false;
        }

        public override Element GetElement()
        {
            Element element = base.GetElement();
            if (element != null)
            {
                return element;
            }
            else if (ContainsConstituentFeatures ())
            {
                //int iCount = ConstituentFeatures().Count;

                //if (AreConstituentFeaturesDirectlyIndexible)
                //    LetConstituentFeaturesBeDirectlyIndexible();

                //for (int i = 0; i < iCount; i++)
                //{
                //    Feature nestedFeature = (Feature)ConstituentFeatures[i];

                //    if ((nestedFeature != null) && (nestedFeature.GetSuccess()))
                //    {
                //        return nestedFeature.Element;
                //    }
                //}
                return null;
            }
            else
            {
                return null;
            }
        }


        public override bool TransformContents
        (
        FeatureUpdateContext updateContext,
        Feature              featureToCopy,
        DTransform3d         transformToUse
        )
        {
            PolygonalGeometry PolygonToCopy = featureToCopy as PolygonalGeometry;

            Type myType = GetType();

            if(PolygonToCopy!=null)
            {
                int iCount = featureToCopy.Count(myType);

                if(iCount>0)
                {

                    DeleteElement(updateContext);  // Indicated by TR 261638.

                    if (featureToCopy.AreConstituentFeaturesDirectlyIndexible)
                        LetConstituentFeaturesBeDirectlyIndexible();

                    for(int i=0; i<iCount; i++)
                    {
                        Feature nestedFeatureToCopy = (Feature)featureToCopy.GetItem((Type) null, i);

                        if((nestedFeatureToCopy!=null)&&(nestedFeatureToCopy.GetSuccess()))
                        {
                            Feature nestedFeature = (Feature)nestedFeatureToCopy.CloneAsChild(this);

                            nestedFeature.AlignOptions(this);

                            if(nestedFeature.TransformContents(updateContext, nestedFeatureToCopy, transformToUse))
                            {
                                nestedFeature.DenoteAsHavingBeenUpdated();

                                nestedFeature.SetSuccess(true);

                                AddConstituentFeatureAndAssignToNamedProperty(updateContext, nestedFeature, nestedFeatureToCopy.Name, i);
                            }
                        }
                    }

                    return true;  // if there are cell contents, then there is no other geometry to transform
                }

                Element elementToCopy = featureToCopy.GetElement();

                if (elementToCopy != null)
                {
                    if (elementToCopy.IsShapeElement())
                    {
                        return base.TransformContents(updateContext, featureToCopy, transformToUse);
                    }
                    else
                    {
                        Element el = DgnTools.CloneElement(elementToCopy);
                        if (el != null)
                        {
                            Transform3d com_transformToUse = DgnTools.ToTransform3d(transformToUse);
                            el.Transform(ref com_transformToUse);
                            //???
                            SetElement(el);
                            return true;
                        }
                    }
                }
            }
            return false;
        }



        /*
        public enum DevelopableOptions
        {
            RuleLines       = MsdDevelopableElementOutputType.msdDevelopableRuleLines,
            PlanarRuleLines = MsdDevelopableElementOutputType.msdDevelopableRuleLinesPlanar,
            Polygons          = MsdDevelopableElementOutputType.msdDevelopablePolygons,
            PlanarPolygons    = MsdDevelopableElementOutputType.msdDevelopablePolygonsPlanar,
            Cones           = MsdDevelopableElementOutputType.msdDevelopableCones,
            PlanarCones     = MsdDevelopableElementOutputType.msdDevelopableConesPlanar
        }
        */

        ///<summary>Approxiate a BSplineSurface with developable panels</summary>
        ///<categories>Surface Based</categories>
        [Update]
        public virtual bool SurfaceDevelopmentPanels
        (
        FeatureUpdateContext updateContext,
        [ParentModel]  CoordinateSystem CoordinateSystem,
        ICurve LeftCurve,
        ICurve RightCurve,
        int SampleCount,
        [Optional, DefaultValue(DevelopableOptions.PlanarPolygons)] DevelopableOptions OutputOption
        )
        {
            DeleteAllSubFeatures(updateContext);
            if((LeftCurve!=null)&&(RightCurve!=null))
            {
                Bentley.Interop.MicroStationDGN.BsplineCurve[] curves = new Bentley.Interop.MicroStationDGN.BsplineCurve[2];

                curves[0] = RightCurve.com_bsplineCurve;
                curves[1] = LeftCurve.com_bsplineCurve;

                SurfaceDevelopmentWithPolygons(updateContext, CoordinateSystem, curves, OutputOption, SampleCount);

                return true;
            }
            return false;
        }
        private bool SurfaceDevelopmentWithPolygons
        (
        FeatureUpdateContext updateContext,
        CoordinateSystem CoordinateSystem,
        Bentley.Interop.MicroStationDGN.BsplineCurve[] curves,
        DevelopableOptions OutputOption,
        int SampleCount
        )
        {
            if ((curves != null) && (curves.Length == 2) && (CoordinateSystem != null))
            {
                /*
                if(  (OutputOption!=DevelopableOptions.PlanarPolygons)
                   &&(OutputOption!=DevelopableOptions.Polygons))
                {
                    OutputOption = DevelopableOptions.PlanarPolygons;
                }
                */

                DTransform3d layoutTransform = CoordinateSystem.m_DTransform3d;

                Transform3d com_layoutTransform = DgnTools.ToTransform3d(layoutTransform);

                ElementEnumerator elementEnumerator = curves[0].ConstructDevelopableRulingsToCurve(curves[1], ref COM_ZERO_POINT3D, (MsdDevelopableElementOutputType)OutputOption, SampleCount);

                if (elementEnumerator != null)
                {
                    Element[] elementArray = (Element[])elementEnumerator.BuildArrayFromContents();

                    if (elementArray != null)
                    {
                        if (elementArray.Length == 1)
                        {
                            if (OutputOption == DevelopableOptions.PlanarPolygons)
                            {
                                com_layoutTransform = DgnTools.ToTransform3d(layoutTransform);

                                elementArray[0].Transform(ref com_layoutTransform);
                            }

                            //???
                            SetElement(elementArray[0]);
                            AddDefaultConstructionBasedOnElement();

                            return true;
                        }

                        else if (elementArray.Length > 1)
                        {
                            int iCount = elementArray.Length;

                            LetConstituentFeaturesBeDirectlyIndexible();

                            for (int i = 0; i < iCount; i++)
                            {
                                Feature nestedFeature;
                                if (OutputOption == DevelopableOptions.RuleLines || OutputOption == DevelopableOptions.PlanarRuleLines)
                                {
                                    nestedFeature = new Line(this);
                                }
                                else
                                {
                                    nestedFeature = new MyFabricationPlanning(this);
                                }

                                nestedFeature.AlignOptions(this);

                                if (OutputOption == DevelopableOptions.PlanarPolygons || OutputOption == DevelopableOptions.PlanarRuleLines)
                                {
                                    com_layoutTransform = DgnTools.ToTransform3d(layoutTransform);

                                    elementArray[i].Transform(ref com_layoutTransform);
                                }

                                //???
                                nestedFeature.SetElement(elementArray[i]);
                                nestedFeature.AddDefaultConstructionBasedOnElement();
                                nestedFeature.DenoteAsHavingBeenUpdated();

                                nestedFeature.SetSuccess(true);

                                SetConstituentFeature(updateContext, i, nestedFeature);
                            }

                            return true;
                        }
                    }
                }
            }

            return false;
        }

        ///<summary>Approxiate a portion of a BSplineSurface with a portion of a cone</summary>
        ///<categories>Surface Based</categories>
        [Update]
        public virtual bool SurfaceDevelopmentCones
        (
        FeatureUpdateContext updateContext,
        [ParentModel] ICurve LeftCurve,
        ICurve RightCurve,
        IPoint SelectorPoint,
        int SampleCount,
        [Out]   ref double[] BaseRadii,
        [Out]   ref double[] TopRadii,
        [Out]   ref Point[] BasePoints,
        [Out]   ref Point[] TopPoints,
        [Out]   ref int selectedIndex
        )
        {
            if((LeftCurve!=null)&&(RightCurve!=null)&&(SelectorPoint!=null))
            {
                DeleteAllSubFeatures(updateContext);
                Bentley.Interop.MicroStationDGN.BsplineCurve[] curves = new Bentley.Interop.MicroStationDGN.BsplineCurve[2];

                curves[0] = RightCurve.com_bsplineCurve;
                curves[1] = LeftCurve.com_bsplineCurve;

                SurfaceDevelopmentWithCones(updateContext,
                                                 curves,
                                                 SelectorPoint,
                                                 DevelopableOptions.Cones,
                                                 SampleCount,
                                                 ref BaseRadii,
                                                 ref TopRadii,
                                                 ref BasePoints,
                                                 ref TopPoints,
                                                 ref selectedIndex);

                return true;
            }

            return false;
        }
        private bool SurfaceDevelopmentWithCones
        (
        FeatureUpdateContext updateContext,
        Bentley.Interop.MicroStationDGN.BsplineCurve[] curves,
        IPoint SelectorPoint,
        DevelopableOptions OutputOption,
        int SampleCount,
        ref double[] BaseRadii,
        ref double[] TopRadii,
        ref Point[] BasePoints,
        ref Point[] TopPoints,
        ref int selectedIndex
        )
        {
            if ((curves != null) && (curves.Length == 2))
            {
                if ((OutputOption != DevelopableOptions.PlanarCones)
                    && (OutputOption != DevelopableOptions.Cones))
                {
                    OutputOption = DevelopableOptions.PlanarCones;
                }

                Point3d origin = DgnTools.ToPoint3d(DPoint3d.Zero);

                ElementEnumerator elementEnumerator = curves[0].ConstructDevelopableRulingsToCurve(curves[1], ref origin, (MsdDevelopableElementOutputType)DevelopableOptions.Polygons, SampleCount);

                if (elementEnumerator != null)
                {

                    Element[] elementArray = (Element[])elementEnumerator.BuildArrayFromContents();
                    MyFabricationPlanning[] polygonalArray = new MyFabricationPlanning[elementArray.Length];
                    double testLength = curves[0].ComputeCurveLength(0, 1, false) + curves[1].ComputeCurveLength(0, 1, false); /* This should, in every case, be longer than the distance between
                                                                                                                            * SelectorPoint and any ruled polygonal */

                    for (int i = 0; i < elementArray.Length; i++)
                    {
                        ShapeElement shapeElement = (ShapeElement)elementArray[i];
                        DPoint3d[] shapeElementVertices = DgnTools.ToDPoint3dArray(shapeElement.GetVertices());
                        DPoint3d pt1 = GeometryTools.CentroidOfPoints(shapeElementVertices);
                        DPoint3d pt2 = SelectorPoint.DPoint3d;

                        if (GeometryTools.Point3dDistance(DgnTools.ToPoint3d(pt1), DgnTools.ToPoint3d(pt2)) < testLength)
                        {
                            testLength = GeometryTools.Point3dDistance(DgnTools.ToPoint3d(pt1), DgnTools.ToPoint3d(pt2));
                            selectedIndex = i;
                        }
                    }

                    //if (ConstructionsVisible == true)
                    //{
                    MyFabricationPlanning nestedFabricationPlanning = new MyFabricationPlanning(this);
                    //???
                    nestedFabricationPlanning.SetElement(elementArray[selectedIndex]);
                    nestedFabricationPlanning.AddDefaultConstructionBasedOnElement();

                    AddConstructionFeature(nestedFabricationPlanning, "RulePolygonTangentToCone");
                    //}
                }

                elementEnumerator = curves[0].ConstructDevelopableRulingsToCurve(curves[1], ref origin, (MsdDevelopableElementOutputType)OutputOption, SampleCount);

                if (elementEnumerator != null)
                {
                    Element[] elementArray = (Element[])elementEnumerator.BuildArrayFromContents();

                    if (elementArray != null)
                    {
                        if (elementArray.Length == 1)
                        {
                            //???
                            SetElement(elementArray[0]);
                            AddDefaultConstructionBasedOnElement();
                        }

                        else
                        {
                            int jCount = elementArray.Length;

                            ConeElement coneElement;
                            TopRadii = new double[jCount];
                            BaseRadii = new double[jCount];
                            TopPoints = new Point[jCount];
                            BasePoints = new Point[jCount];
                            for (int j = 0; j < jCount; j++)
                            {
                                coneElement = elementArray[j].AsConeElement();

                                TopPoints[j] = new Point(this);
                                TopPoints[j].FromDPoint3d(DgnTools.ToDPoint3d(coneElement.get_TopCenterPoint()));
                                BasePoints[j] = new Point(this);
                                BasePoints[j].FromDPoint3d(DgnTools.ToDPoint3d(coneElement.get_BaseCenterPoint()));
                                TopRadii[j] = Math.Abs(coneElement.TopRadius);
                                BaseRadii[j] = Math.Abs(coneElement.BaseRadius);

                            }

                            //???
                            SetElement(elementArray[selectedIndex]);
                            AddDefaultConstructionBasedOnElement();

                            UpdateDisplay();

                            /*
                            for(int i=iCount-1; i>-1; i--)
                            {
                                SetElement(DgnTools.CloneElement(elementArray[i]););

                                UpdateDisplay();

                                Thread.Sleep(TimeInterval);
                            }
                            */

                        }
                    }
                }

                return true;
            }

            return false;
        }

    }
}
