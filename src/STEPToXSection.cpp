//
// Program: STEPToXSection
//
// Description:
//
// The program STEPToXSection extracts the contour of a planar cross section of solids contained in STEP files.
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
#include <BRepAlgoAPI_Section.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <gp_Pln.hxx>
#include <GCPnts_UniformDeflection.hxx>
#include <vector>
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

auto computeXSection(const TopoDS_Compound& compound, const std::array<double, 4>& plane, const double deflection) -> std::vector<std::array<double, 3>> {
	std::vector<std::array<double, 3>> points;
	BRepAlgoAPI_Section section(compound, gp_Pln(plane[0], plane[1], plane[2], plane[3]), Standard_False);
	section.ComputePCurveOn1(Standard_True);
	section.Approximation(Standard_True);
	section.Build();
	std::list<TopoDS_Edge> edges;
	TopExp_Explorer explorer;
	for (explorer.Init(section.Shape(), TopAbs_EDGE); explorer.More(); explorer.Next()) {
		TopoDS_Edge edge{TopoDS::Edge(explorer.Current())};
		BRepAdaptor_Curve curve{edge};
		GCPnts_UniformDeflection discretizer(curve, deflection, Standard_True);
		if (discretizer.IsDone() && discretizer.NbPoints() > 0) {
			int nbPoints = discretizer.NbPoints();
			for (int i = 2; i <= nbPoints; i++) {
				const gp_Pnt p0 = discretizer.Value(i - 1);
				const gp_Pnt p1 = discretizer.Value(i);
				points.emplace_back(std::array<double, 3>{p0.X(), p0.Y(), p0.Z()});
				points.emplace_back(std::array<double, 3>{p1.X(), p1.Y(), p1.Z()});
			}
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

void write(const std::string& outFile, const std::vector<NamedSolid>& namedSolids, const std::vector<std::string>& names,
		const double deflection, const std::array<double, 4>& plane, const std::string& format) {
	TopoDS_Compound compound;
	TopoDS_Builder builder;
	builder.MakeCompound(compound);
	if (names.empty()) {
		for (const auto& namedSolid : namedSolids) builder.Add(compound, namedSolid.solid);
	}
	else {
		for (const auto& name : names) {
			const auto iter = std::find_if(std::begin(namedSolids), std::end(namedSolids), [&](const auto& namesSolid) { return namesSolid.name == name; });
			if (iter == std::end(namedSolids)) throw std::logic_error{std::string{"Could not find solid with name '"} + name + "'"};
			builder.Add(compound, iter->solid);
		}
	}
	std::vector<std::array<double, 3>> points{computeXSection(compound, plane, deflection)};
	if (format == "xyz") writeXYZ(outFile, plane, points);
	else if (format == "ply") writePLY(outFile, plane, points);
}

int main(int argc, char* argv[]) {
	cxxopts::Options options{"STEPToXSection", "Extracts the contour of a planar cross section of solids contained in the STEP file"};
	options.add_options()
			("i,in", "Input file", cxxopts::value<std::string>())
			("o,out", "Output file", cxxopts::value<std::string>())
			("f,format", "Output file format (xyz or ply)", cxxopts::value<std::string>()->default_value("xyz"))
			("c,content", "List content (solids)")
			("s,select", "Select solids by name (comma seperated list)", cxxopts::value<std::vector<std::string>>())
			("d,deflection", "deflection", cxxopts::value<double>())
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
			if (plane.size() != 4) throw std::logic_error{std::string{"Wrong plane format (nx,ny,nz,d)"}};
			const std::array<double, 4> planeVec{plane[0], plane[1], plane[2], plane[3]};
			const auto format{result["format"].as<std::string>()};
			if (format != "xyz" && format != "ply") throw std::logic_error{std::string{"Format '"} + format + "' not supported"};
			const auto deflection{result["deflection"].as<double>()};
			std::vector<std::string> select;
			if (result.count("select")) select = result["select"].as<std::vector<std::string>>();
			std::vector<NamedSolid> namedSolids;
			read(inFile, namedSolids);
			write(outFile, namedSolids, select, deflection, planeVec, format);
		}
		else std::cout << options.help() << std::endl;
		return EXIT_SUCCESS;
	}
	catch (const std::exception& ex) {
		std::cout << ex.what();
		return EXIT_FAILURE;
	}
}
