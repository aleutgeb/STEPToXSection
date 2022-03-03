//
// Program: STEPToXSection
//
// Description:
//
// The program STEPToXSection extracts the contour of a planar cross section of solids contained in STEP files. Additionally it supports offsetting of the cross section contour.
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
#include <TopExp_Explorer.hxx>
#include <STEPCAFControl_Reader.hxx>
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <TDocStd_Document.hxx>
#include <TDF_ChildIterator.hxx>
#include <TDataStd_Name.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
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
#include <BRep_Tool.hxx>
#include <BRepFeat.hxx>
#include <GeomLProp_SLProps.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <gp_Pln.hxx>
#include <GCPnts_UniformDeflection.hxx>
#include <vector>
#include <unordered_map>
#include <deque>
#include <array>
#include <cmath>
#include <iostream>
#include "cxxopts.hpp"

struct NamedSolid {
	NamedSolid(const TopoDS_Solid& s, const std::string& n) : solid{s}, name{n} {}

	const TopoDS_Solid solid;
	const std::string  name;
};

void getNamedSolids(const TopLoc_Location& location, const std::string& prefix, unsigned int& id, const Handle(XCAFDoc_ShapeTool) shapeTool,
		const TDF_Label label, std::vector<NamedSolid>& namedSolids) {
	TDF_Label referredLabel{label};
	if (shapeTool->IsReference(label)) shapeTool->GetReferredShape(label, referredLabel);
	std::string name;
	Handle(TDataStd_Name) shapeName;
	if (referredLabel.FindAttribute(TDataStd_Name::GetID(), shapeName)) name = TCollection_AsciiString(shapeName->Get()).ToCString();
	if (name == "") name = std::to_string(id++);
	std::string fullName{prefix + "/" + name};

	TopLoc_Location localLocation = location * shapeTool->GetLocation(label);
	TDF_LabelSequence components;
	if (shapeTool->GetComponents(referredLabel, components)) {
		for (Standard_Integer compIndex = 1; compIndex <= components.Length(); ++compIndex) {
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
	Handle(XCAFApp_Application) application = XCAFApp_Application::GetApplication();
	application->NewDocument(inFile.c_str(), document);
	STEPCAFControl_Reader reader;
	reader.SetNameMode(true);
	IFSelect_ReturnStatus stat = reader.ReadFile(inFile.c_str());
	if (stat != IFSelect_RetDone || !reader.Transfer(document)) throw std::logic_error{std::string{"Could not read '"} + inFile + "'"};
	Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(document->Main());
	TDF_LabelSequence topLevelShapes;
	shapeTool->GetFreeShapes(topLevelShapes);
	unsigned int id{1};
	for (Standard_Integer iLabel = 1; iLabel <= topLevelShapes.Length(); ++iLabel) {
		getNamedSolids(TopLoc_Location{}, "", id, shapeTool, topLevelShapes.Value(iLabel), namedSolids);
	}
}

auto convertToPolygon(const TopoDS_Wire& wire, const gp_Pln& xsectionPlane, const double deflection) -> TopoDS_Wire {
	BRepBuilderAPI_MakePolygon makePolygon;
	BRepTools_WireExplorer wireExplorer;
	for (wireExplorer.Init(wire); wireExplorer.More(); wireExplorer.Next()) {
		TopoDS_Edge edge{wireExplorer.Current()};
		BRepAdaptor_Curve curve{edge};
		GCPnts_UniformDeflection discretizer{curve, deflection, Standard_True};
		if (discretizer.IsDone() && discretizer.NbPoints() > 0) {
			int nbPoints = discretizer.NbPoints();
			std::vector<gp_Pnt> points;
			for (int i{1}; i <= nbPoints; i++) {
				gp_Pnt p{discretizer.Value(i)};
				// p lies on xsection plane, but BRepOffsetAPI_MakeOffset needs a more exact positioning. So we move p closer to the plane.
				gp_XYZ coord{p.Coord()};
				coord += xsectionPlane.Axis().Direction().XYZ().Multiplied(xsectionPlane.Distance(p));
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
	for (explorer.Init(shape, TopAbs_EDGE); explorer.More(); explorer.Next())
	{
		edges->Append(TopoDS::Edge(explorer.Current()));
	}
	ShapeAnalysis_FreeBounds::ConnectEdgesToWires(edges, 0.0, Standard_True, wires);
	for (int iWire{1}; iWire <= wires->Size(); ++iWire) result.emplace_back(TopoDS::Wire(wires->Value(iWire)));
	return result;
}

auto computeXSectionWires(const TopoDS_Compound& compound, const gp_Pln& xsectionPlane) -> std::vector<TopoDS_Wire> {
	BRepAlgoAPI_Section section(compound, xsectionPlane, Standard_False);
	section.ComputePCurveOn1(Standard_True);
	section.Approximation(Standard_True);
	section.Build();
	return convertToWires(section.Shape());
}

void discretizeWire(const TopoDS_Wire& wire, std::vector<std::array<double, 3>>& points, const double deflection) {
	TopExp_Explorer explorer;
	for (explorer.Init(wire, TopAbs_EDGE); explorer.More(); explorer.Next()) {
		TopoDS_Edge edge{TopoDS::Edge(explorer.Current())};
		BRepAdaptor_Curve curve{edge};
		GCPnts_UniformDeflection discretizer(curve, deflection, Standard_True);
		if (discretizer.IsDone() && discretizer.NbPoints() > 0) {
			int nbPoints = discretizer.NbPoints();
			for (int i = 2; i <= nbPoints; i++) {
				const gp_Pnt p0{discretizer.Value(i - 1)};
				const gp_Pnt p1{discretizer.Value(i)};
				points.emplace_back(std::array<double, 3>{p0.X(), p0.Y(), p0.Z()});
				points.emplace_back(std::array<double, 3>{p1.X(), p1.Y(), p1.Z()});
			}
		}
	}
}

auto surfaceNormal(const TopoDS_Face& face) -> gp_Dir {
	Standard_Real umin, umax, vmin, vmax;
	BRepTools::UVBounds(face, umin, umax, vmin, vmax);
	Handle(Geom_Surface) surface {BRep_Tool::Surface(face)};
	GeomLProp_SLProps props(surface, umin, vmin, 1, 0.01);
	return props.Normal();
}

auto cutTwoFaces(const TopoDS_Face& origFace0, const TopoDS_Face& origFace1) -> std::vector<TopoDS_Face> {
	std::vector<TopoDS_Face> result;
	TopTools_ListOfShape arguments;
	TopTools_ListOfShape tools;
	TopoDS_Face face0{origFace0};
	TopoDS_Face face1{origFace1};
	if (surfaceNormal(face0).Dot(surfaceNormal(face1)) < 0.0) face1.Reverse();
	arguments.Append(face0);
	tools.Append(face1);
	BRepAlgoAPI_Cut cut;
	cut.SetArguments(arguments);
	cut.SetTools(tools);
	cut.Build();
	TopoDS_Shape shape = cut.Shape();
	ShapeUpgrade_UnifySameDomain unify(shape, Standard_True, Standard_True, Standard_False);
	unify.Build();
	shape = unify.Shape();
	for (TopoDS_Iterator iter(shape); iter.More(); iter.Next()) {
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
	TopoDS_Shape shape = fuse.Shape();
	ShapeUpgrade_UnifySameDomain unify(shape, Standard_True, Standard_True, Standard_False);
	unify.Build();
	shape = unify.Shape();
	for (TopoDS_Iterator iter(shape); iter.More(); iter.Next()) {
		if (iter.Value().ShapeType() == TopAbs_FACE) {
			result.emplace_back(TopoDS::Face(iter.Value()));
		}
	}
	return result;
}

void uniteFaces(const std::vector<TopoDS_Face>& faces, std::vector<std::array<double, 3>>& points, const double deflection) {
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
		if (!united) {
			for (const auto& wire : convertToWires(face)) discretizeWire(wire, points, deflection);
		}
	}
}

auto computeContainment(const std::vector<TopoDS_Wire>& wires) -> std::unordered_map<int, std::vector<int>> {
	std::unordered_map<int, std::vector<int>> result;
	std::vector<bool> contained(wires.size(), false);
	std::vector<std::vector<int>> contains(wires.size());
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

auto computeOffsetWire(const TopoDS_Wire& wire, const double offset) -> std::optional<TopoDS_Wire> {
	std::optional<TopoDS_Wire> result;
	BRepOffsetAPI_MakeOffset makeOffset{wire, GeomAbs_Arc};
	makeOffset.Perform(offset);
	if (makeOffset.IsDone()) {
		const TopoDS_Shape& shape{makeOffset.Shape()};
		if (!shape.IsNull() && shape.ShapeType() == TopAbs_WIRE) {
			result = TopoDS::Wire(shape);
		}
	}
	else throw std::runtime_error{ "Invalid offset operation" };
	return result;
}

auto computeXSection(const TopoDS_Compound& compound, const std::array<double, 4>& plane, const double deflection, const double offset) -> std::vector<std::array<double, 3>> {
	std::vector<std::array<double, 3>> points;
	const gp_Pln xsectionPlane{plane[0], plane[1], plane[2], plane[3]};
	std::vector<TopoDS_Wire> wires{computeXSectionWires(compound, xsectionPlane)};
	if (std::abs(offset) < 0.1 * deflection) {
		for (const auto& wire : wires) discretizeWire(wire, points, deflection);
	}
	else {
		std::vector<TopoDS_Wire> polygonWires;
		for (const auto& wire : wires) polygonWires.emplace_back(convertToPolygon(wire, xsectionPlane, deflection));
		std::unordered_map<int, std::vector<int>> containment{computeContainment(polygonWires)};
		std::vector<TopoDS_Face> offsetFaces;
		for (const auto entry : containment) {
			std::optional<TopoDS_Wire> outerWire{computeOffsetWire(polygonWires[entry.first], offset)};
			if (outerWire) {
				BRepBuilderAPI_MakeFace makeFace{*outerWire};
				TopoDS_Face face{makeFace.Face()};
				for (const auto inner : entry.second) {
					std::optional<TopoDS_Wire> innerWire{computeOffsetWire(polygonWires[inner], -offset) };
					if (innerWire) {
						std::vector<TopoDS_Face> cut{cutTwoFaces(face, BRepBuilderAPI_MakeFace{*innerWire}.Face())};
						if (cut.size() == 1) face = cut[0];
					}
				}
				offsetFaces.emplace_back(face);
			}
		}
		if (offsetFaces.size() == 1) {
			for (const auto& wire : convertToWires(offsetFaces[0])) discretizeWire(wire, points, deflection);
		}
		else if (offsetFaces.size() > 1) {
			uniteFaces(offsetFaces, points, deflection);
		}
	}
	return points;
}

void writeXYZ(const std::string& outFile, const std::array<double, 4>& plane, const std::vector<std::array<double, 3>>& points) {
	std::ofstream ofs{outFile};
	for (const auto& p : points) {
		ofs << p[0] << " " << p[1] << " " << p[2] << " " << plane[0] << " " << plane[1] << " " << plane[2] << "\n";
	}
	ofs.close();
}

void writePLY(const std::string& outFile, const std::array<double, 4>& plane, const std::vector<std::array<double, 3>>& points) {
	std::ofstream ofs{outFile};
	ofs << "ply" << "\n";
	ofs << "format ascii 1.0" << "\n";
	ofs << "element vertex " << points.size() << "\n";
	ofs << "property float x" << "\n";
	ofs << "property float y" << "\n";
	ofs << "property float z" << "\n";
	ofs << "element edge " << (points.size() / 2) << "\n";
	ofs << "property int vertex1" << "\n";
	ofs << "property int vertex2" << "\n";
	ofs << "end_header" << "\n";
	for (std::size_t i{0}; i < points.size(); ++i) {
		const auto& p{points[i]};
		ofs << p[0] << " " << p[1] << " " << p[2] << "\n";
	}
	for (std::size_t i{0}; i < points.size(); i += 2) {
		ofs << i << " " << i + 1 << "\n";
	}
	ofs.close();
}

void write(const std::string& outFile, const std::vector<NamedSolid>& namedSolids, const std::vector<std::string>& select,
		const double deflection, const double offset, const std::array<double, 4>& plane, const std::string& format) {
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
					const auto iter = std::find_if(std::begin(namedSolids), std::end(namedSolids), [&](const auto& namesSolid) { return namesSolid.name == sel; });
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
	std::vector<std::array<double, 3>> points{computeXSection(compound, plane, deflection, offset)};
	if (format == "xyz") writeXYZ(outFile, plane, points);
	else if (format == "ply") writePLY(outFile, plane, points);
}

int main(int argc, char* argv[]) {
	cxxopts::Options options{"STEPToXSection", "Extracts the contour of a planar cross section of solids contained in the STEP file. Additionally it supports offsetting of the cross section contour."};
	options.add_options()
			("i,in", "Input file", cxxopts::value<std::string>())
			("o,out", "Output file", cxxopts::value<std::string>())
			("f,format", "Output file format (xyz or ply)", cxxopts::value<std::string>()->default_value("xyz"))
			("c,content", "List content (solids)")
			("s,select", "Select solids by name or index (comma seperated list, index starts with 1)", cxxopts::value<std::vector<std::string>>())
			("d,deflection", "deflection", cxxopts::value<double>())
			("t,offset", "offset", cxxopts::value<double>()->default_value("0.0"))
			("p,plane", "Plane (a,b,c,d), in which a*x + b*y + c*z + d = 0", cxxopts::value<std::vector<double>>())
			("h,help", "Print usage");
	try {
		const auto result = options.parse(argc, argv);
		if (result.count("content")) {
			if (result.count("in")) {
				const std::string inFile = result["in"].as<std::string>();
				std::vector<NamedSolid> namedSolids;
				read(inFile, namedSolids);
				for (const auto& namedSolid : namedSolids) std::cout << namedSolid.name << std::endl;
			}
			else throw std::logic_error{std::string{"Missing option 'in'"}};
		}
		else if (result.count("in") && result.count("out")) {
			const auto inFile{result["in"].as<std::string>()}, outFile{result["out"].as<std::string>()};
			if (!result.count("deflection")) throw std::logic_error{std::string{"Missing option 'deflection'"}};
			if (!result.count("plane")) throw std::logic_error{std::string{"Missing option 'plane'"}};
			std::vector<double> plane{result["plane"].as<std::vector<double>>()};
			if (plane.size() != 4) throw std::logic_error{std::string{"Wrong plane format (a,b,c,d)"}};
			const std::array<double, 4> planeVec{plane[0], plane[1], plane[2], plane[3]};
			const auto format{result["format"].as<std::string>()};
			if (format != "xyz" && format != "ply") throw std::logic_error{std::string{"Format '"} + format + "' not supported"};
			const auto deflection{result["deflection"].as<double>()};
			const auto offset{result["offset"].as<double>()};
			std::vector<std::string> select;
			if (result.count("select")) select = result["select"].as<std::vector<std::string>>();
			std::vector<NamedSolid> namedSolids;
			read(inFile, namedSolids);
			write(outFile, namedSolids, select, deflection, offset, planeVec, format);
		}
		else std::cout << options.help() << std::endl;
		return EXIT_SUCCESS;
	}
	catch (const std::exception& ex) {
		std::cout << ex.what();
		return EXIT_FAILURE;
	}
}
