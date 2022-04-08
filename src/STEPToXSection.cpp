//
// Program: STEPToXSection
//
// Description:
//
// The program STEPToXSection extracts the contour of a planar cross section of solids contained in STEP files.
// It supports the orthogonal projection of geometries with a specified maximum plane distance, in which the silhouette of the projected geometries represents the base contour.
// Additionally the program supports the computation of contour offset curves.
//
// Copyright(C) 2021 Alexander Leutgeb
//
// This library is free software; you can redistribute it and / or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110 - 1301  USA
//

#include <TopoDS_Solid.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Builder.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopExp_Explorer.hxx>
#include <STEPCAFControl_Reader.hxx>
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <TDocStd_Document.hxx>
#include <TDF_ChildIterator.hxx>
#include <TDataStd_Name.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <ShapeAnalysis_FreeBounds.hxx>
#include <BRepAlgoAPI_Section.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepOffsetAPI_MakeOffset.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRepPrimAPI_MakeHalfSpace.hxx>
#include <BRep_Tool.hxx>
#include <BRepFeat.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <GeomLProp_SLProps.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <ShapeFix_Wire.hxx>
#include <gp_Pln.hxx>
#include <gp_Trsf.hxx>
#include <HLRBRep_Algo.hxx>
#include <HLRBRep_HLRToShape.hxx>
#include <GCPnts_UniformDeflection.hxx>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <deque>
#include <numeric>
#include <array>
#include <cmath>
#include <iostream>
#include "cxxopts.hpp"

struct NamedSolid {
	NamedSolid(const TopoDS_Solid& s, const std::string& n) : solid{s}, name{n} {}

	const TopoDS_Solid solid;
	const std::string  name;
};

struct ResultEntry {
	unsigned int plane;
	unsigned int offset;
	TopoDS_Wire  wire;
};

struct pair_hash {
	template <class T1, class T2>
	std::size_t operator() (const std::pair<T1, T2>& pair) const {
		return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
	}
};

struct Connectivity {
	std::vector<gp_XYZ>                                         vertices;
	std::unordered_map<unsigned int, std::vector<unsigned int>> vertexToVertices;
};

struct CoordinatesLessOperator {
	CoordinatesLessOperator(const double e) : eps{e} {}

	bool operator()(const gp_XYZ& lhs, const gp_XYZ& rhs) const {
		for (auto i{0u}; i < 3u; ++i) {
			if (std::abs(lhs.GetData()[i] - rhs.GetData()[i]) > eps) return lhs.GetData()[i] < rhs.GetData()[i];
		}
		return false;
	}

	const double eps;
};

void writeXYZ(const std::string& outFile, const std::array<double, 3>& planeNormal, const std::vector<ResultEntry>& result) {
	std::ofstream ofs{outFile};
	for (const auto& entry : result) {
		TopExp_Explorer explorer;
		for (explorer.Init(entry.wire, TopAbs_EDGE); explorer.More(); explorer.Next()) {
			const TopoDS_Edge edge{TopoDS::Edge(explorer.Current())};
			const gp_XYZ p0{BRep_Tool::Pnt(TopExp::FirstVertex(edge)).Coord()};
			const gp_XYZ p1{BRep_Tool::Pnt(TopExp::LastVertex(edge)).Coord()};
			ofs << p0.X() << " " << p0.Y() << " " << p0.Z() << " " << planeNormal[0] << " " << planeNormal[1] << " " << planeNormal[2] << "\n";
			ofs << p1.X() << " " << p1.Y() << " " << p1.Z() << " " << planeNormal[0] << " " << planeNormal[1] << " " << planeNormal[2] << "\n";
		}
	}
	ofs.close();
}

void writePLY(const std::string& outFile, const std::vector<TopoDS_Edge>& edges) {
	std::ofstream ofs{outFile};
	ofs << "ply" << "\n";
	ofs << "format ascii 1.0" << "\n";
	ofs << "element vertex " << edges.size() * 2 << "\n";
	ofs << "property double x" << "\n";
	ofs << "property double y" << "\n";
	ofs << "property double z" << "\n";
	ofs << "element edge " << edges.size() << "\n";
	ofs << "property int vertex1" << "\n";
	ofs << "property int vertex2" << "\n";
	ofs << "end_header" << "\n";
	for (const TopoDS_Edge& edge : edges) {
		const std::array<gp_XYZ, 2> points{BRep_Tool::Pnt(TopExp::FirstVertex(edge)).Coord(), BRep_Tool::Pnt(TopExp::LastVertex(edge)).Coord()};
		for (const auto& p : points) ofs << p.X() << " " << p.Y() << " " << p.Z() << "\n";
	}
	for (std::size_t i{0}; i < edges.size(); ++i) ofs << i * 2 << " " << i * 2 + 1 << "\n";
	ofs.close();
}

void writePLYEdges(const std::string& outFile, const std::vector<ResultEntry>& result) {
	std::vector<TopoDS_Edge> edges;
	for (const auto& entry : result) {
		TopExp_Explorer explorer;
		for (explorer.Init(entry.wire, TopAbs_EDGE); explorer.More(); explorer.Next()) edges.emplace_back(TopoDS::Edge(explorer.Current()));
	}
	writePLY(outFile, edges);
}

void writePLY(const std::string& outFile, const std::vector<ResultEntry>& result) {
	std::vector<std::array<double, 3>> vertices;
	std::vector<std::vector<unsigned int>> faces;
	for (const auto& entry : result) {
		faces.push_back({});
		TopExp_Explorer explorer;
		for (explorer.Init(entry.wire, TopAbs_EDGE); explorer.More(); explorer.Next()) {
			const TopoDS_Edge edge{TopoDS::Edge(explorer.Current())};
			const gp_XYZ p{BRep_Tool::Pnt(TopExp::FirstVertex(edge)).Coord()};
			faces.back().emplace_back(static_cast<unsigned int>(vertices.size()));
			vertices.emplace_back(std::array<double, 3>{p.X(), p.Y(), p.Z()});
		}
	}

	std::ofstream ofs{outFile};
	ofs << "ply" << "\n";
	ofs << "format ascii 1.0" << "\n";
	ofs << "element vertex " << vertices.size() << "\n";
	ofs << "property double x" << "\n";
	ofs << "property double y" << "\n";
	ofs << "property double z" << "\n";
	ofs << "element face " << faces.size() << "\n";
	ofs << "property list int int vertex_index" << "\n";
	ofs << "property int plane_index" << "\n";
	ofs << "property int offset_index" << "\n";
	ofs << "end_header" << "\n";
	for (const auto& v : vertices) ofs << v[0] << " " << v[1] << " " << v[2] << "\n";
	for (std::size_t i{0u}; i < faces.size(); ++i) {
		const auto& face{faces[i]};
		ofs << face.size();
		for (const auto& vIdx : face) ofs << " " << vIdx;
		ofs << " " << result[i].plane << " " << result[i].offset << "\n";
	}
	ofs.close();
}

void getNamedSolids(const TopLoc_Location& location, const std::string& prefix, unsigned int& id, const Handle(XCAFDoc_ShapeTool) shapeTool,
		const TDF_Label label, std::vector<NamedSolid>& namedSolids) {
	TDF_Label referredLabel{label};
	if (shapeTool->IsReference(label)) shapeTool->GetReferredShape(label, referredLabel);
	std::string name;
	Handle(TDataStd_Name) shapeName;
	if (referredLabel.FindAttribute(TDataStd_Name::GetID(), shapeName)) name = TCollection_AsciiString(shapeName->Get()).ToCString();
	if (name == "") name = std::to_string(id++);
	std::string fullName{prefix + "/" + name};

	TopLoc_Location localLocation{location * shapeTool->GetLocation(label)};
	TDF_LabelSequence components;
	if (shapeTool->GetComponents(referredLabel, components)) {
		for (Standard_Integer compIndex{1}; compIndex <= components.Length(); ++compIndex) {
			getNamedSolids(localLocation, fullName, id, shapeTool, components.Value(compIndex), namedSolids);
		}
	}
	else {
		TopoDS_Shape shape;
		shapeTool->GetShape(referredLabel, shape);
		if (shape.ShapeType() == TopAbs_SOLID) {
			BRepBuilderAPI_Transform transform(shape, localLocation, Standard_True);
			namedSolids.emplace_back(TopoDS::Solid(transform.Shape()), fullName);
		}
	}
}

void read(const std::string& inFile, std::vector<NamedSolid>& namedSolids) {
	Handle(TDocStd_Document) document;
	Handle(XCAFApp_Application) application{XCAFApp_Application::GetApplication()};
	application->NewDocument(inFile.c_str(), document);
	STEPCAFControl_Reader reader;
	reader.SetNameMode(true);
	IFSelect_ReturnStatus stat{reader.ReadFile(inFile.c_str())};
	if (stat != IFSelect_RetDone || !reader.Transfer(document)) throw std::logic_error{std::string{"Could not read '"} + inFile + "'"};
	Handle(XCAFDoc_ShapeTool) shapeTool{XCAFDoc_DocumentTool::ShapeTool(document->Main())};
	TDF_LabelSequence topLevelShapes;
	shapeTool->GetFreeShapes(topLevelShapes);
	unsigned int id{1u};
	for (Standard_Integer iLabel{1}; iLabel <= topLevelShapes.Length(); ++iLabel) {
		getNamedSolids(TopLoc_Location{}, "", id, shapeTool, topLevelShapes.Value(iLabel), namedSolids);
	}
}

auto convertToPolygonOnPlane(const TopoDS_Wire& wire, const gp_Pln& xsectionPlane, const double deflection) -> TopoDS_Wire {
	BRepBuilderAPI_MakePolygon makePolygon;
	BRepTools_WireExplorer wireExplorer;
	for (wireExplorer.Init(wire); wireExplorer.More(); wireExplorer.Next()) {
		TopoDS_Edge edge{wireExplorer.Current()};
		BRepAdaptor_Curve curve{edge};
		GCPnts_UniformDeflection discretizer{curve, deflection, Standard_True};
		if (discretizer.IsDone() && discretizer.NbPoints() > 0) {
			int nbPoints{discretizer.NbPoints()};
			std::vector<gp_Pnt> points;
			for (int i{1}; i <= nbPoints; i++) {
				gp_Pnt p{discretizer.Value(i)};
				gp_XYZ coord{p.Coord()};
				const gp_XYZ vec{p.XYZ() - xsectionPlane.Location().XYZ()};
				const double distance{vec * xsectionPlane.Axis().Direction().XYZ()};
				coord -= xsectionPlane.Axis().Direction().XYZ().Multiplied(distance);
				p.SetCoord(coord.X(), coord.Y(), coord.Z());
				points.emplace_back(p.X(), p.Y(), p.Z());
			}
			const gp_Pnt startPoint{BRep_Tool::Pnt(wireExplorer.CurrentVertex())};
			if (points.size() > 0) {
				if (points.back().Distance(startPoint) < points.front().Distance(startPoint)) std::reverse(std::begin(points), std::end(points));
				for (auto i{1u}; i < points.size(); ++i) makePolygon.Add(points[i]);
			}
		}
	}
	makePolygon.Close();
	return makePolygon.Wire();
}

auto convertToWires(const TopoDS_Shape& shape) -> std::vector<TopoDS_Wire> {
	std::vector<TopoDS_Wire> result;
	Handle(TopTools_HSequenceOfShape) edges{new TopTools_HSequenceOfShape{}};
	Handle(TopTools_HSequenceOfShape) wires{new TopTools_HSequenceOfShape{}};
	TopExp_Explorer explorer;
	for (explorer.Init(shape, TopAbs_EDGE); explorer.More(); explorer.Next()) edges->Append(TopoDS::Edge(explorer.Current()));
	ShapeAnalysis_FreeBounds::ConnectEdgesToWires(edges, 0.001, Standard_False, wires);
	for (int iWire{1}; iWire <= wires->Size(); ++iWire) result.emplace_back(TopoDS::Wire(wires->Value(iWire)));
	return result;
}

auto computeXSectionWires(const TopoDS_Compound& compound, const gp_Pln& xsectionPlane) -> std::vector<TopoDS_Wire> {
	BRepAlgoAPI_Section section{compound, xsectionPlane, Standard_False};
	section.ComputePCurveOn2(Standard_True);
	section.Approximation(Standard_True);
	section.Build();
	return convertToWires(section.Shape());
}

auto surfaceNormal(const TopoDS_Face& face) -> gp_Dir {
	Standard_Real umin, umax, vmin, vmax;
	BRepTools::UVBounds(face, umin, umax, vmin, vmax);
	Handle(Geom_Surface) surface{BRep_Tool::Surface(face)};
	GeomLProp_SLProps props{surface, umin, vmin, 1, 0.01};
	return props.Normal();
}

auto cutFaces(const TopoDS_Face& face, const std::vector<TopoDS_Face>& faces) -> std::vector<TopoDS_Face> {
	std::vector<TopoDS_Face> result;
	TopTools_ListOfShape arguments;
	arguments.Append(face);
	TopTools_ListOfShape tools;
	for (const auto& f : faces) {
		TopoDS_Face newF{f};
		if (surfaceNormal(face).Dot(surfaceNormal(newF)) < 0.0) newF.Reverse();
		tools.Append(newF);
	}
	BRepAlgoAPI_Cut cut;
	cut.SetArguments(arguments);
	cut.SetTools(tools);
	cut.Build();
	TopoDS_Shape shape{cut.Shape()};
	ShapeUpgrade_UnifySameDomain unify{shape, Standard_True, Standard_True, Standard_False};
	unify.Build();
	shape = unify.Shape();
	for (TopoDS_Iterator iter{shape}; iter.More(); iter.Next()) {
		if (iter.Value().ShapeType() == TopAbs_FACE) {
			result.emplace_back(TopoDS::Face(iter.Value()));
		}
	}
	return result;
}

auto uniteTwoFaces(const TopoDS_Face& origFace0, const TopoDS_Face& origFace1) -> std::vector<TopoDS_Face> {
	std::vector<TopoDS_Face> result;
	TopTools_ListOfShape arguments;
	TopTools_ListOfShape tools;
	TopoDS_Face face0{origFace0};
	TopoDS_Face face1{origFace1};
	if (surfaceNormal(face0).Dot(surfaceNormal(face1)) < 0.0) face1.Reverse();
	arguments.Append(face0);
	tools.Append(face1);
	BRepAlgoAPI_Fuse fuse;
	fuse.SetArguments(arguments);
	fuse.SetTools(tools);
	fuse.Build();
	TopoDS_Shape shape{fuse.Shape()};
	ShapeUpgrade_UnifySameDomain unify{shape, Standard_True, Standard_True, Standard_False};
	unify.Build();
	shape = unify.Shape();
	for (TopoDS_Iterator iter{shape}; iter.More(); iter.Next()) {
		if (iter.Value().ShapeType() == TopAbs_FACE) {
			result.emplace_back(TopoDS::Face(iter.Value()));
		}
	}
	return result;
}

auto uniteFaces(const std::vector<TopoDS_Face>& faces) -> std::vector<TopoDS_Face> {
	std::vector<TopoDS_Face> result;
	std::deque<TopoDS_Face> work;
	std::copy(std::begin(faces), std::end(faces), std::back_inserter(work));
	while (!work.empty()) {
		TopoDS_Face face{work.front()};
		work.pop_front();
		bool united{false};
		for (auto iter{std::begin(work)}; iter != std::end(work); ++iter) {
			std::vector<TopoDS_Face> result{uniteTwoFaces(face, *iter)};
			if (result.size() == 1) {
				united = true;
				work.erase(iter);
				work.push_back(result[0]);
				break;
			}
		}
		if (!united) result.emplace_back(face);
	}
	return result;
}

auto computeContainment(const std::vector<TopoDS_Wire>& wires) -> std::unordered_map<unsigned int, std::vector<unsigned int>> {
	std::unordered_map<unsigned int, std::vector<unsigned int>> result;
	std::vector<bool> contained(wires.size(), false);
	for (auto iOuter{0u}; iOuter < wires.size(); ++iOuter) {
		TopoDS_Face outerFace{BRepBuilderAPI_MakeFace{wires[iOuter]}.Face()};
		for (auto iInner{0u}; iInner < wires.size(); ++iInner) {
			if (iOuter != iInner) {
				TopoDS_Face innerFace{BRepBuilderAPI_MakeFace{wires[iInner]}.Face()};
				if (BRepFeat::IsInside(innerFace, outerFace)) {
					contained[iInner] = true;
					result[iOuter].emplace_back(iInner);
				}
			}
		}
	}
	for (auto iOuter{0u}; iOuter < wires.size(); ++iOuter) {
		if (!contained[iOuter] && result.find(iOuter) == std::end(result)) result[iOuter];
	}
	return result;
}

auto splitAtFirstSelfIntersection(const TopoDS_Wire& wire) -> std::vector<TopoDS_Wire> {
	std::vector<TopoDS_Wire> result;
	std::vector<TopoDS_Edge> edges;
	std::vector<gp_Pnt> points;

	BRepTools_WireExplorer wireExplorer;
	for (wireExplorer.Init(wire); wireExplorer.More(); wireExplorer.Next()) {
		edges.emplace_back(wireExplorer.Current());
		points.emplace_back(BRep_Tool::Pnt(wireExplorer.CurrentVertex()));
	}

	ShapeAnalysis_Wire analysis{wire, BRepBuilderAPI_MakeFace{wire}.Face(), 0.0};
	std::optional<std::tuple<unsigned int, unsigned int, gp_Pnt>> intersection;
	for (auto i{1u}; i + 1u <= edges.size() && !intersection; ++i) {
		for (auto j{i + 1u}; j <= edges.size() && !intersection; ++j) {
			if (analysis.CheckIntersectingEdges(i, j)) {
				IntRes2d_SequenceOfIntersectionPoint points2d;
				TColgp_SequenceOfPnt points3d;
				TColStd_SequenceOfReal errors;
				analysis.CheckIntersectingEdges(i, j, points2d, points3d, errors);
				if (points3d.Size() > 0) intersection = std::make_tuple(i, j, points3d.Value(1));
			}
		}
	}
	if (intersection) {
		auto [start, end, point]{*intersection};
		BRepBuilderAPI_MakePolygon makePolygon0;
		for (auto i{start + 1u}; i <= end; ++i) makePolygon0.Add(points[i - 1]);
		makePolygon0.Add(point);
		makePolygon0.Close();
		result.emplace_back(makePolygon0.Wire());
		BRepBuilderAPI_MakePolygon makePolygon1;
		for (auto i{end + 1u}; i <= points.size() + start; ++i) makePolygon1.Add(points[(i - 1) % points.size()]);
		makePolygon1.Add(point);
		makePolygon1.Close();
		result.emplace_back(makePolygon1.Wire());
	}
	else result.emplace_back(wire);
	return result;
}

auto splitWireAtSelfIntersections(const TopoDS_Wire& wire, const TopoDS_Face& face) -> std::vector<TopoDS_Wire> {
	std::vector<TopoDS_Wire> result;
	std::vector<TopoDS_Wire> allWires;
	std::deque<TopoDS_Wire> work;
	work.emplace_back(wire);
	while (!work.empty()) {
		const auto currentWire{work.front()};
		work.pop_front();
		std::vector<TopoDS_Wire> subWires{splitAtFirstSelfIntersection(currentWire)};
		if (subWires.size() == 1) allWires.emplace_back(subWires[0]);
		else {
			for (const auto& subWire : subWires) work.emplace_back(subWire);
		}
	}
	for (const auto& w : allWires) {
		if (surfaceNormal(face).Dot(surfaceNormal(BRepBuilderAPI_MakeFace{w}.Face())) > 0.0) result.emplace_back(w);
	}
	return result;
}

auto computeOffsetWires(const TopoDS_Wire& wire, const double offset) -> std::vector<TopoDS_Wire> {
	std::vector<TopoDS_Wire> result;
	BRepOffsetAPI_MakeOffset makeOffset{wire, GeomAbs_Arc};
	makeOffset.Perform(offset);
	if (makeOffset.IsDone()) {
		const TopoDS_Shape shape{makeOffset.Shape()};
		if (!shape.IsNull()){
			if (shape.ShapeType() == TopAbs_WIRE) {
				if (offset > 0.0) result.emplace_back(TopoDS::Wire(TopoDS::Wire(shape).Reversed()));
				else result.emplace_back(TopoDS::Wire(shape));
			}
		}
		else throw std::runtime_error{"Wire offset computation failed"};
	}
	else throw std::runtime_error{"Wire offset computation failed"};
	return result;
}

auto computeOffsetPolygonWires(const TopoDS_Wire& wire, const double offset, const gp_Pln& xsectionPlane, const double deflection) -> std::vector<TopoDS_Wire> {
	std::vector<TopoDS_Wire> result{computeOffsetWires(wire, offset)};
	for (auto i{0u}; i < result.size(); ++i) result[i] = convertToPolygonOnPlane(result[i], xsectionPlane, deflection);
	return result;
}

auto computeSplittedOffsetPolygonWires(const TopoDS_Wire& wire, const double offset, const gp_Pln& xsectionPlane, const double deflection) -> std::vector<TopoDS_Wire> {
	std::vector<TopoDS_Wire> result;
	if (std::abs(offset) < 0.1 * deflection) result.emplace_back(wire);
	else {
		std::vector<TopoDS_Wire> offsetWires{computeOffsetPolygonWires(wire, offset, xsectionPlane, deflection)};
		for (const auto& offsetWire : offsetWires) {
			for (const auto& splittedWire : splitWireAtSelfIntersections(offsetWire, BRepBuilderAPI_MakeFace{wire}.Face())) result.emplace_back(splittedWire);
		}
	}
	return result;
}

auto discretizeEdge(const TopoDS_Edge& edge, const gp_Trsf& transformation, const gp_Pln& xsectionPlane, double deflection) -> std::vector<TopoDS_Edge> {
	std::vector<TopoDS_Edge> result;
	BRepAdaptor_Curve curve{edge};
	GCPnts_UniformDeflection discretizer{curve, deflection, Standard_True};
	if (discretizer.IsDone() && discretizer.NbPoints() > 0) {
		int nbPoints{discretizer.NbPoints()};
		for (int i{2}; i <= nbPoints; i++) {
			std::array<gp_Pnt, 2> points{discretizer.Value(i - 1), discretizer.Value(i)};
			for (auto& p : points) {
				p.Transform(transformation);
				gp_XYZ coord{p.Coord()};
				const gp_XYZ vec{p.XYZ() - xsectionPlane.Location().XYZ()};
				const double distance{vec * xsectionPlane.Axis().Direction().XYZ()};
				coord -= xsectionPlane.Axis().Direction().XYZ().Multiplied(distance);
				p.SetCoord(coord.X(), coord.Y(), coord.Z());
			}
			result.emplace_back(BRepBuilderAPI_MakeEdge(points[0], points[1]).Edge());
		}
	}
	return result;
}

auto projectShape(const TopoDS_Shape& shape, const gp_Pln& xsectionPlane, const double distance, const double deflection) -> std::vector<TopoDS_Edge> {
	gp_Pln plane0{xsectionPlane.Location().XYZ() - xsectionPlane.Axis().Direction().XYZ() * distance, -xsectionPlane.Axis().Direction()};
	gp_Pln plane1{xsectionPlane.Location().XYZ() + xsectionPlane.Axis().Direction().XYZ() * distance, xsectionPlane.Axis().Direction()};
	const TopoDS_Face face0{BRepBuilderAPI_MakeFace{plane0}};
	const TopoDS_Face face1{BRepBuilderAPI_MakeFace{plane1}};
	const TopoDS_Shape halfSpace0{BRepPrimAPI_MakeHalfSpace{face0, xsectionPlane.Location().XYZ() - xsectionPlane.Axis().Direction().XYZ() * distance * 2.0}.Solid()};
	const TopoDS_Shape halfSpace1{BRepPrimAPI_MakeHalfSpace{face1, xsectionPlane.Location().XYZ() + xsectionPlane.Axis().Direction().XYZ() * distance * 2.0}.Solid()};
	TopTools_ListOfShape arguments;
	arguments.Append(shape);
	TopTools_ListOfShape tools;
	tools.Append(halfSpace0);
	tools.Append(halfSpace1);
	BRepAlgoAPI_Cut cut;
	cut.SetFuzzyValue(deflection);
	cut.SetArguments(arguments);
	cut.SetTools(tools);
	cut.Build();
	TopoDS_Shape cutShape{cut.Shape()};
	const gp_Ax2 axis{xsectionPlane.Location().XYZ() - xsectionPlane.Axis().Direction().XYZ() * distance * 2.0, xsectionPlane.Axis().Direction()};
	Handle(HLRBRep_Algo) hlr{new HLRBRep_Algo};
	hlr->Add(cutShape);
	HLRAlgo_Projector projector{axis};
	const gp_Trsf inverse{projector.InvertedTransformation()};
	hlr->Projector(projector);
	hlr->Update();
	hlr->Hide();
	HLRBRep_HLRToShape hlrToShape(hlr);
	std::vector<TopoDS_Edge> result;
	for (const auto type : { /*HLRBRep_Undefined, HLRBRep_IsoLine,*/ HLRBRep_OutLine, HLRBRep_Rg1Line, /*HLRBRep_RgNLine , */HLRBRep_Sharp}) {
		TopoDS_Shape hlrShape{hlrToShape.CompoundOfEdges(type, true, false)};
		for (TopoDS_Iterator iter{hlrShape}; iter.More(); iter.Next()) {
			if (iter.Value().ShapeType() == TopAbs_EDGE) {
				const std::vector<TopoDS_Edge> edges{discretizeEdge(TopoDS::Edge(iter.Value()), inverse, xsectionPlane, deflection)};
				std::copy(std::begin(edges), std::end(edges), std::back_inserter(result));
			}
		}
	}
	return result;
}

auto computeXSectionFaces(const std::vector<TopoDS_Wire>& wires, const gp_Pln& xsectionPlane, const double deflection, const double offset) -> std::vector<TopoDS_Face> {
	std::vector<TopoDS_Face> result;
	std::vector<TopoDS_Wire> polygonWires;
	for (const auto& wire : wires) polygonWires.emplace_back(convertToPolygonOnPlane(wire, xsectionPlane, deflection));
	std::unordered_map<unsigned int, std::vector<unsigned int>> containment{computeContainment(polygonWires)};
	std::vector<TopoDS_Face> offsetFaces;
	for (const auto& entry : containment) {
		std::vector<TopoDS_Wire> outerWires{computeSplittedOffsetPolygonWires(polygonWires[entry.first], offset, xsectionPlane, deflection)};
		for (const auto& oW : outerWires) {
			TopoDS_Face face{BRepBuilderAPI_MakeFace{oW}.Face()};
			std::vector<TopoDS_Face> cuttingFaces;
			for (const auto inner : entry.second) {
				std::vector<TopoDS_Wire> innerWires{computeSplittedOffsetPolygonWires(polygonWires[inner], -offset, xsectionPlane, deflection)};
				for (const auto& iW : innerWires) cuttingFaces.emplace_back(BRepBuilderAPI_MakeFace{iW}.Face());
			}
			if (cuttingFaces.empty()) offsetFaces.emplace_back(face);
			else for (const auto& f : cutFaces(face, cuttingFaces)) offsetFaces.emplace_back(f);
		}
	}
	return uniteFaces(offsetFaces);
}

void disconnectEdge(Connectivity& connectivity, const unsigned int from, const unsigned int to) {
	auto& forwards{connectivity.vertexToVertices[from]};
	forwards.erase(std::find(std::begin(forwards), std::end(forwards), to));
	auto& backwards{connectivity.vertexToVertices[to]};
	backwards.erase(std::find(std::begin(backwards), std::end(backwards), from));
}

void connectToVertices(Connectivity& connectivity, const double eps) {
	std::vector<std::pair<unsigned int, unsigned int>> joinEdges;
	for (const auto& fromTo : connectivity.vertexToVertices) {
		if (fromTo.second.size() < 2) {
			const auto& fromPoint{connectivity.vertices[fromTo.first]};
			auto minDist{std::numeric_limits<double>::max()};
			unsigned int minIdx;
			for (auto pointIdx{0u}; pointIdx < connectivity.vertices.size(); ++pointIdx) {
				if (pointIdx != fromTo.first) {
					const auto dist{(connectivity.vertices[pointIdx] - fromPoint).Modulus()};
					if (dist < minDist) {
						minDist = dist;
						minIdx = pointIdx;
					}
				}
			}
			if (minDist < eps) joinEdges.emplace_back(fromTo.first, minIdx);
		}
	}
	for (const auto& fromTo : joinEdges) {
		connectivity.vertexToVertices[fromTo.first].emplace_back(fromTo.second);
		connectivity.vertexToVertices[fromTo.second].emplace_back(fromTo.first);
	}
}

void connectToEdges(Connectivity& connectivity, const double eps) {
	std::unordered_map<std::pair<unsigned int, unsigned int>, std::vector<unsigned int>, pair_hash> splitEdges;
	for (const auto& fromTo : connectivity.vertexToVertices) {
		if (fromTo.second.size() < 2) {
			const auto& fromPoint{connectivity.vertices[fromTo.first]};
			auto minDist{std::numeric_limits<double>::max()};
			unsigned int minFromIdx;
			unsigned int minToIdx;
			for (const auto& vertexEdges : connectivity.vertexToVertices) {
				const auto& edgeFrom{vertexEdges.first};
				const auto& edgeTos{vertexEdges.second};
				if (fromTo.first != edgeFrom) {
					for (const auto& edgeTo : edgeTos) {
						if (fromTo.first != edgeTo) {
							BRepBuilderAPI_MakeEdge makeEdge{connectivity.vertices[edgeFrom], connectivity.vertices[edgeTo]};
							if (makeEdge.IsDone()) {
								BRepExtrema_DistShapeShape distCalc;
								distCalc.LoadS1(makeEdge.Edge());
								distCalc.LoadS2(BRepBuilderAPI_MakeVertex{fromPoint}.Shape());
								distCalc.Perform();
								const auto dist{distCalc.Value()};
								if (dist < minDist) {
									minDist = dist;
									minFromIdx = edgeFrom;
									minToIdx = edgeTo;
								}
							}
						}
					}
				}
			}
			if (minDist < eps) {
				auto edge{std::make_pair(minFromIdx, minToIdx)};
				if (edge.first > edge.second) std::swap(edge.first, edge.second);
				splitEdges[edge].emplace_back(fromTo.first);
			}
		}
	}

	for (const auto& edgeToPoints : splitEdges) {
		const auto& from{edgeToPoints.first.first};
		const auto& to{edgeToPoints.first.second};
		const auto& point{edgeToPoints.second};
		std::vector<unsigned int> sorted{edgeToPoints.second};
		std::sort(std::begin(sorted), std::end(sorted), [&connectivity, &from](const auto idx0, const auto idx1) {
			const auto& p{connectivity.vertices[from]};
			const auto& p0{connectivity.vertices[idx0]};
			const auto& p1{connectivity.vertices[idx1]};
			return (p0 - p).Modulus() < (p1 - p).Modulus();
		});
		disconnectEdge(connectivity, from, to);
		sorted.insert(std::begin(sorted), from);
		sorted.emplace_back(to);
		for (auto i{0u}; i + 1u < sorted.size(); ++i) {
			connectivity.vertexToVertices[sorted[i]].emplace_back(sorted[i + 1u]);
			connectivity.vertexToVertices[sorted[i + 1u]].emplace_back(sorted[i]);
		}
	}
}

void connectOpenEdges(Connectivity& connectivity, const double eps) {
	connectToVertices(connectivity, eps * 5.0);
	connectToEdges(connectivity, eps);
}

void writeConnectivity(const Connectivity& connectivity, const std::string& outFile) {
	std::vector<TopoDS_Edge> edges;
	for (const auto& fromTos : connectivity.vertexToVertices) {
		const auto& p0{connectivity.vertices[fromTos.first]};
		for (const auto& to : fromTos.second) {
			const auto& p1{connectivity.vertices[to]};
			BRepBuilderAPI_MakeEdge makeEdge{p0, p1};
			if (makeEdge.IsDone()) edges.emplace_back(makeEdge.Edge());
		}
	}
	writePLY(outFile, edges);
}

auto computeConnectivity(const std::vector<TopoDS_Edge>& origEdges, const double eps) -> Connectivity {
	Connectivity result;
	std::vector<std::pair<unsigned int, unsigned int>> edges;

	CoordinatesLessOperator comp{eps};
	std::map <gp_XYZ, unsigned int, CoordinatesLessOperator> vertexMap{comp};
	for (auto iEdge{0u}; iEdge < origEdges.size(); ++iEdge) {
		const auto& edge{origEdges[iEdge]};
		const std::array<gp_XYZ, 2> points{BRep_Tool::Pnt(TopExp::FirstVertex(edge)).Coord(), BRep_Tool::Pnt(TopExp::LastVertex(edge)).Coord()};
		std::array<unsigned int, 2> indices;
		for (auto iSide{0u}; iSide < 2u; ++iSide) {
			auto it{vertexMap.find(points[iSide])};
			if (it == std::end(vertexMap)) {
				indices[iSide] = static_cast<unsigned int>(result.vertices.size());
				vertexMap[points[iSide]] = indices[iSide];
				result.vertices.emplace_back(points[iSide]);
			}
			else indices[iSide] = it->second;
		}
		edges.emplace_back(indices[0], indices[1]);
	}

	for (const auto& fromTo : edges) {
		result.vertexToVertices[fromTo.first].emplace_back(fromTo.second);
		result.vertexToVertices[fromTo.second].emplace_back(fromTo.first);
	}
	connectOpenEdges(result, eps);
	return result;
}

void computeCluster(const std::unordered_map<unsigned int, std::vector<unsigned int>>& vertexToVertices,
		const unsigned int current, std::unordered_set<unsigned int>& visited, std::vector<unsigned int>& cluster) {
	cluster.emplace_back(current);
	visited.insert(current);
	auto it{vertexToVertices.find(current)};
	if (it != std::end(vertexToVertices)) {
		for (const auto& next : it->second) {
			if (visited.find(next) == std::end(visited)) {
				computeCluster(vertexToVertices, next, visited, cluster);
			}
		}
	}
}

auto findSeed(const Connectivity& connectivity, const std::vector<unsigned int>& cluster, const double deflection) -> std::optional<unsigned int> {
	std::array<double, 3> min{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
	std::array<double, 3> max{std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};

	for (const auto& idx : cluster) {
		const gp_XYZ& p{connectivity.vertices[idx]};
		for (auto dim{0u}; dim < 3u; ++dim) {
			min[dim] = std::min(min[dim], p.GetData()[dim]);
			max[dim] = std::max(max[dim], p.GetData()[dim]);
		}
	}

	std::array<unsigned int, 3> dims{0u, 1u, 2u};
	std::sort(std::begin(dims), std::end(dims), [&min, &max](const auto lhs, const auto rhs) { return max[lhs] - min[lhs] > max[rhs] - min[rhs]; });
	std::vector<unsigned int> candidates{cluster};
	for (auto dim{0u}; dim < 3u && candidates.size() > 1; ++dim) {
		const auto curDim{dims[dim]};
		std::vector<unsigned int> newCandidates;
		std::sort(std::begin(candidates), std::end(candidates), [&connectivity, &curDim](const auto idx0, const auto idx1) {
			const gp_XYZ& p0{connectivity.vertices[idx0]};
			const gp_XYZ& p1{connectivity.vertices[idx1]};
			return p0.GetData()[curDim] < p1.GetData()[curDim];
		});
		const gp_XYZ& ref{connectivity.vertices[candidates.front()]};
		for (const auto& idx : candidates) {
			const gp_XYZ& p{connectivity.vertices[idx]};
			if (std::abs(p.GetData()[curDim] - ref.GetData()[curDim]) > deflection) break;
			else newCandidates.emplace_back(idx);
		}
		candidates = newCandidates;
	}
	if (!candidates.empty()) return candidates.front();
	else return std::nullopt;
}

auto findCandidate(const Connectivity& connectivity, const gp_XYZ& normal, const unsigned int prev, const unsigned int current) -> std::optional<unsigned int> {
	const auto c_minAngle{std::acos(-1.0) * 1.0e-6};
	const std::vector<unsigned int>& nexts{connectivity.vertexToVertices.find(current)->second};
	const auto& prevPoint{connectivity.vertices[prev]};
	const auto& currentPoint{connectivity.vertices[current]};
	const auto edgePrev{(prevPoint - currentPoint).Normalized()};
	double maxAngle{0.0};
	std::optional<unsigned int> result;
	for (const auto& next : nexts) {
		if (next != prev) {
			const auto& nextPoint{connectivity.vertices[next]};
			const auto edgeNext{(nextPoint - currentPoint).Normalized()};
			auto angle{std::acos(std::clamp(edgePrev * edgeNext, -1.0, 1.0))};
			if (angle > c_minAngle && angle + c_minAngle < std::acos(-1.0) && (edgeNext ^ edgePrev).Normalized() * normal < 0.0) angle = 2.0 * std::acos(-1.0) - angle;
			if (angle > maxAngle) {
				maxAngle = angle;
				result = next;
			}
		}
	}
	return result;
}

auto computeWire(const Connectivity& connectivity, const std::vector<unsigned int>& cluster, const double deflection) -> std::optional<TopoDS_Wire> {
	std::optional<unsigned int> currentOpt{findSeed(connectivity, cluster, deflection)};
	if (!currentOpt) return std::nullopt;
	unsigned int current{*currentOpt};
	auto currentPoint{connectivity.vertices[current]};
	auto tos{connectivity.vertexToVertices.find(current)->second};
	auto minDP{std::numeric_limits<double>::max()};
	unsigned int prev;
	unsigned int next;
	for (auto i{0u}; i < tos.size(); ++i) {
		const auto& toPoint0{connectivity.vertices[tos[i]]};
		const auto edge0{(toPoint0 - currentPoint).Normalized()};
		for (auto j{i + 1u}; j < tos.size(); ++j) {
			const auto& toPoint1{connectivity.vertices[tos[j]]};
			const auto edge1{(toPoint1 - currentPoint).Normalized()};
			const auto dp{edge0 * edge1};
			if (dp < minDP) {
				minDP = dp;
				next = tos[i];
				prev = tos[j];
			}
		}
	}
	if (!(minDP < 1.0)) return std::nullopt;
	std::vector<unsigned int> loop;
	std::unordered_set<unsigned int> visited;
	loop.emplace_back(prev);
	loop.emplace_back(current);
	loop.emplace_back(next);
	visited.insert(prev);
	visited.insert(current);
	visited.insert(next);
	const auto& prevPoint{connectivity.vertices[prev]};
	const auto& nextPoint{connectivity.vertices[next]};
	const auto edgeNext{(nextPoint - currentPoint).Normalized()};
	const auto edgePrev{(prevPoint - currentPoint).Normalized()};
	const auto normal{(edgeNext ^ edgePrev).Normalized()};
	while (true) {
		auto candidate{findCandidate(connectivity, normal, current, next)};
		if (!candidate || visited.find(*candidate) != std::end(visited)) break;
		loop.emplace_back(*candidate);
		visited.insert(*candidate);
		current = next;
		next = *candidate;
	}
	BRepBuilderAPI_MakePolygon makePolygon;
	for (const auto idx : loop) makePolygon.Add(connectivity.vertices[idx]);
	makePolygon.Close();
	return makePolygon.Wire();
}

auto computeWires(const Connectivity& connectivity, const double deflection) -> std::vector<TopoDS_Wire> {
	std::vector<TopoDS_Wire> result;
	std::unordered_set<unsigned int> visited;
	std::vector<unsigned int> vertexIndices;
	for (const auto& fromTos : connectivity.vertexToVertices) vertexIndices.emplace_back(fromTos.first);
	auto it{std::find_if(std::begin(vertexIndices), std::end(vertexIndices),[&visited](const unsigned int index) { return visited.find(index) == std::end(visited); })};
	while (it != std::end(vertexIndices)) {
		std::vector<unsigned int> cluster;
		computeCluster(connectivity.vertexToVertices, *it, visited, cluster);
		const auto wire{computeWire(connectivity, cluster, deflection)};
		if (wire) result.emplace_back(*wire);
		it = std::find_if(std::begin(vertexIndices), std::end(vertexIndices), [&visited](const unsigned int index) { return visited.find(index) == std::end(visited); });
	}
	return result;
}

auto computeProjectedSilhouette(const TopoDS_Compound& compound, const gp_Pln& xsectionPlane, const double projection, const double deflection) -> std::vector<TopoDS_Wire> {
	const std::vector<TopoDS_Edge> edges{projectShape(compound, xsectionPlane, projection, deflection)};
	Connectivity connectivity{computeConnectivity(edges, deflection)};
	return computeWires(connectivity, deflection);
}

auto computeXSection(const TopoDS_Compound& compound, const std::array<double, 3>& planeNormal, const double planeDistance, const double deflection,
		const double offset, const double projection) -> std::vector<TopoDS_Wire> {
	std::vector<TopoDS_Wire> result;
	const gp_Pln xsectionPlane{planeNormal[0], planeNormal[1], planeNormal[2], planeDistance};
	std::vector<TopoDS_Wire> wires;
	if (projection > deflection * 0.1) wires = computeProjectedSilhouette(compound, xsectionPlane, projection, deflection);
	else wires = computeXSectionWires(compound, xsectionPlane);
	std::vector<TopoDS_Face> faces{computeXSectionFaces(wires, xsectionPlane, deflection, offset)};
	for (const auto& face : faces) {
		for (const auto& wire : convertToWires(face)) result.emplace_back(wire);
	}
	return result;
}

auto computeXSections(const TopoDS_Compound& compound, const std::array<double, 3>& planeNormal, const std::vector<double>& planeDistances,
		const double deflection, const std::vector<double>& offsets, const double projection) -> std::vector<ResultEntry> {
	std::vector<ResultEntry> result;
	for (auto iPlane{0u}; iPlane < planeDistances.size(); ++iPlane) {
		const auto& planeDistance{planeDistances[iPlane]};
		for (auto iOffset{0u}; iOffset < offsets.size(); ++iOffset) {
			const auto& offset{offsets[iOffset]};
			const auto wires{computeXSection(compound, planeNormal, planeDistance, deflection, offset, projection)};
			for (const auto& wire : wires) result.emplace_back(ResultEntry{iPlane, iOffset, wire});
		}
	}
	return result;
}

auto getSelectedSolids(const std::vector<NamedSolid>& namedSolids, const std::vector<std::string>& select) -> TopoDS_Compound {
	TopoDS_Compound compound;
	TopoDS_Builder builder;
	builder.MakeCompound(compound);
	if (select.empty()) {
		for (const auto& namedSolid : namedSolids) builder.Add(compound, namedSolid.solid);
	}
	else {
		for (const auto& sel : select) {
			if (sel != "") {
				if (sel[0] == '/') {
					const auto iter{std::find_if(std::begin(namedSolids), std::end(namedSolids), [&](const auto& namesSolid) { return namesSolid.name == sel; })};
					if (iter == std::end(namedSolids)) throw std::logic_error{std::string{"Could not find solid with name '"} + sel + "'"};
					builder.Add(compound, iter->solid);
				}
				else {
					try {
						int index{std::stoi(sel)};
						if (index < 1 || index > namedSolids.size()) throw std::logic_error{std::string{"Index out of range: "} + sel};
						builder.Add(compound, namedSolids[index - 1].solid);
					}
					catch (const std::invalid_argument&) {
						throw std::logic_error{std::string("Invalid index: ") + sel};
					}
				}
			}
		}
	}
	return compound;
}

void write(const std::string& outFile, const std::vector<NamedSolid>& namedSolids, const std::vector<std::string>& select,
		const double deflection, const std::vector<double>& offsets, const std::array<double, 3>& planeNormal, const std::vector<double>& planeDistances,
		const double projection, const std::string& format) {
	TopoDS_Compound compound{getSelectedSolids(namedSolids, select)};
	const std::vector<ResultEntry> result{computeXSections(compound, planeNormal, planeDistances, deflection, offsets, projection)};
	if (format == "xyz") writeXYZ(outFile, planeNormal, result);
	else if (format == "ply") writePLY(outFile, result);
}

auto generateRange(const double start, const double end, const double count) -> std::vector<double> {
	std::vector<double> result;
	const auto iCount{static_cast<int>(std::floor(count))};
	if (iCount < 2) throw std::logic_error{std::string{"Count must be greater 1"}};
	const auto delta{end - start};
	const auto inc{delta / (iCount - 1)};
	auto value{start};
	for (auto i{0}; i < iCount; ++i, value += inc) result.emplace_back(value);
	return result;
}

int main(int argc, char* argv[]) {
	cxxopts::Options options{ "STEPToXSection", "Extracts the contour of a planar cross section of solids contained in STEP files. "
			"It supports the orthogonal projection of geometries with a specified maximum plane distance, in which the silhouette of the projected geometries represents the base contour. "
			"Additionally the program supports the computation of contour offset curves."};
	options.add_options()
			("i,in", "Input file", cxxopts::value<std::string>())
			("o,out", "Output file", cxxopts::value<std::string>())
			("f,format", "Output file format (xyz or ply)", cxxopts::value<std::string>()->default_value("ply"))
			("c,content", "List content (solids)")
			("s,select", "Select solids by name or index (comma seperated list, index starts with 1)", cxxopts::value<std::vector<std::string>>())
			("d,deflection", "Chordal tolerance used during discretization", cxxopts::value<double>())
			("p,plane", "Single plane (a,b,c,d) or parallel planes (a,b,c,d_start,d_end,d_count), in which a*x + b*y + c*z + d = 0", cxxopts::value<std::vector<double>>())
			("t,offset", "Single offset (value) or range offset (start,end,count) for computation of contour offset curves", cxxopts::value<std::vector<double>>())
			("n,projection", "Orthogonal projection of geometries with specified maximum plane distance, in which the silhouette of the projected geometries represents the base contour", cxxopts::value<double>())
			("h,help", "Print usage");
	try {
		const auto result{options.parse(argc, argv)};
		if (result.count("content")) {
			if (result.count("in")) {
				const std::string inFile{result["in"].as<std::string>()};
				std::vector<NamedSolid> namedSolids;
				read(inFile, namedSolids);
				for (const auto& namedSolid : namedSolids) std::cout << namedSolid.name << std::endl;
			}
			else throw std::logic_error{std::string{"Missing option 'in'"}};
		}
		else if (result.count("in") && result.count("out")) {
			const auto inFile{result["in"].as<std::string>()}, outFile{result["out"].as<std::string>()};
			const auto format{result["format"].as<std::string>()};
			if (format != "xyz" && format != "ply") throw std::logic_error{std::string{"Format '"} + format + "' not supported"};
			std::vector<std::string> select;
			if (result.count("select")) select = result["select"].as<std::vector<std::string>>();
			if (!result.count("deflection")) throw std::logic_error{std::string{"Missing option 'deflection'"}};
			const auto deflection{result["deflection"].as<double>()};
			if (!result.count("plane")) throw std::logic_error{std::string{"Missing option 'plane'"}};
			const std::vector<double> planeParams{result["plane"].as<std::vector<double>>()};
			if (planeParams.size() != 4 && planeParams.size() != 6) throw std::logic_error{std::string{"Wrong plane format, must be (a,b,c,d) or (a,b,c,d_start,d_end,d_count)"}};
			std::vector<double> planeDistances;
			const std::array<double, 3> planeNormal{planeParams[0], planeParams[1], planeParams[2]};
			if (planeParams.size() == 4) planeDistances.emplace_back(planeParams[3]);
			else planeDistances = generateRange(planeParams[3], planeParams[4], planeParams[5]);
			std::vector<double> offsets;
			if (!result.count("offset")) offsets.emplace_back(0.0);
			else {
				std::vector<double> offsetParams{result["offset"].as<std::vector<double>>()};
				if (offsetParams.size() != 1 && offsetParams.size() != 3) throw std::logic_error{std::string{"Invalid offset, must be single value or range (start,end,count)"}};
				if (offsetParams.size() == 1) offsets.emplace_back(offsetParams[0]);
				else offsets = generateRange(offsetParams[0], offsetParams[1], offsetParams[2]);
			}
			double projection{0.0};
			if (result.count("projection")) {
				projection = result["projection"].as<double>();
				if (projection < 0.0) throw std::logic_error{std::string{"Invalid projection, must be >= 0.0"}};
			}
			std::vector<NamedSolid> namedSolids;
			read(inFile, namedSolids);
			write(outFile, namedSolids, select, deflection, offsets, planeNormal, planeDistances, projection, format);
		}
		else std::cout << options.help() << std::endl;
		return EXIT_SUCCESS;
	}
	catch (const std::exception& ex) {
		std::cerr << ex.what() << std::endl;;
		return EXIT_FAILURE;
	}
	catch (...) {
		std::cerr << "Unexpected exception" << std::endl;
		return EXIT_FAILURE;
	}
}
