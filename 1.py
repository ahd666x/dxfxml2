# این اسکریپت برای خواندن و پردازش فایل‌های DXF به کتابخانه 'ezdxf' نیاز دارد.
# لطفاً قبل از اجرا، این کتابخانه را نصب کنید (مثلاً با استفاده از: pip install ezdxf)

import ezdxf
import xml.etree.ElementTree as ET
import xml.dom.minidom # Import minidom for pretty printing
import math
import numpy as np
import re
import os
import sys
import subprocess # Import subprocess module
import unittest # Import the unittest framework

# --- Configuration ---
# Define DXF layer names used by the script
# This makes the script easier to update if layer names change
DXF_LAYER_CONFIG = {
    'part_border': '_ABF_PART_BORDER',
    'cutting_lines': '_ABF_CUTTING_LINES', # Still needed for front-only panels
    'sheet_border': '_ABF_SHEET_BORDER', # Layer for sheet borders
    'drilling_prefix': 'ABF_D', # Prefix for drilling layers (e.g., ABF_D10, ABF_D15)
    'pocket_layer': 'ABF_DSIDE_8', # Layer for pocket/side holes
    'groove_layer': 'ABF_GROOVE' # Layer for grooves
}

# --- Helper Functions ---

def get_bbox_center(vertices):
    """Calculates the center point of a bounding box from vertices."""
    if not vertices: return None, None
    min_x = min(v[0] for v in vertices)
    max_x = max(v[0] for v in vertices)
    min_y = min(v[1] for v in vertices)
    max_y = max(v[1] for v in vertices)
    center_x = (min_x + max_x) / 2.0
    center_y = (min_y + max_y) / 2.0
    return center_x, center_y

def get_bbox(vertices):
    """Calculates the bounding box coordinates from vertices."""
    min_x = float('inf')
    max_x = -float('inf')
    min_y_val = float('inf')
    max_y_val = -float('inf')

    if not vertices: return min_x, min_y_val, max_x, max_y_val # Handle empty vertices

    for v in vertices:
        min_x = min(min_x, v[0])
        max_x = max(max_x, v[0])
        min_y_val = min(min_y_val, v[1])
        max_y_val = max(max_y_val, v[1])

    return min_x, min_y_val, max_x, max_y_val


def get_bbox_dimensions_sorted(vertices):
    """Calculates bounding box dimensions and returns them sorted (smallest, largest)."""
    min_x, min_y, max_x, max_y = get_bbox(vertices)
    if min_x == float('inf'): return 0.0, 0.0 # Handle empty vertices
    width = max_x - min_x
    height = max_y - min_y
    return round(min(width, height), 3), round(max(width, height), 3)

# Removed unused distance function
# def distance(p1, p2):
#     """Calculates the Euclidean distance between two 2D points."""
#     if p1 is None or p2 is None: return float('inf')
#     return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)

def is_point_inside_bbox(point, bbox, tolerance=0.1):
    """Checks if a point is inside a bounding box with tolerance."""
    if point is None or bbox is None: return False
    return (bbox[0] - tolerance <= point[0] <= bbox[2] + tolerance and
            bbox[1] - tolerance <= point[1] <= bbox[3] + tolerance)


def get_number_from_layer_name_after_D(layer_name, drilling_prefix):
    """Extracts the number after the specified drilling prefix from layer names."""
    # Ensure the prefix is treated as literal text in the regex
    escaped_prefix = re.escape(drilling_prefix)
    match = re.search(rf'{escaped_prefix}(\d+)', layer_name.upper())
    if match:
        try:
            return int(match.group(1))
        except ValueError:
            return None
    return None

# --- Core Logic Functions ---

def find_and_group_panels(doc, config):
    """
    Finds _ABF_SHEET_BORDER, _ABF_PART_BORDER, and _ABF_CUTTING_LINES entities
    and groups them into physical panels based on containment within sheet borders
    and the presence of _ABF_PART_BORDER. Now only matches panels based on size,
    not position symmetry.
    Uses layer names from config.
    """
    sheet_border_layer = config['sheet_border'].upper()
    part_border_layer = config['part_border'].upper()
    cutting_lines_layer = config['cutting_lines'].upper()

    sheet_borders = [e for e in doc.modelspace() if e.dxftype() == 'LWPOLYLINE' and e.dxf.layer.upper() == sheet_border_layer]
    part_borders = [e for e in doc.modelspace() if e.dxftype() == 'LWPOLYLINE' and e.dxf.layer.upper() == part_border_layer]
    cutting_lines_borders = [e for e in doc.modelspace() if e.dxftype() == 'LWPOLYLINE' and e.dxf.layer.upper() == cutting_lines_layer]

    if len(sheet_borders) != 2:
        print(f"❌ Error: انتظار می شود دقیقاً دو موجودیت در لایه '{sheet_border_layer}' وجود داشته باشد، اما {len(sheet_borders)} یافت شد. نمی توان پنل ها را گروه بندی کرد.")
        return []

    # Determine which sheet border is "back" (contains part borders) and which is "front"
    sheet_border_bboxes = [get_bbox(list(sb.vertices())) for sb in sheet_borders]
    front_sheet_border_entity = None
    back_sheet_border_entity = None
    front_sheet_bbox = None
    back_sheet_bbox = None

    # Find the sheet border that contains part borders
    for i, sb_bbox in enumerate(sheet_border_bboxes):
        contains_part_border = False
        for pb in part_borders:
            pb_center = get_bbox_center(list(pb.vertices()))
            if pb_center and is_point_inside_bbox(pb_center, sb_bbox, tolerance=1.1):
                contains_part_border = True
                break
        if contains_part_border:
            back_sheet_border_entity = sheet_borders[i]
            back_sheet_bbox = sb_bbox
        else:
            front_sheet_border_entity = sheet_borders[i]
            front_sheet_bbox = sb_bbox

    if not front_sheet_border_entity or not back_sheet_border_entity:
         print(f"❌ Error: نتوانست مرزهای ورق جلو و پشت را بر اساس محتوای لایه '{part_border_layer}' شناسایی کند.")
         return []

    print(f"DEBUG: مرز ورق جلو شناسایی شد (محدوده: {front_sheet_bbox})")
    print(f"DEBUG: مرز ورق پشت شناسایی شد (محدوده: {back_sheet_bbox})")

    grouped_panels = []
    used_cutting_lines_indices = set()

    # 1. Process panels with _ABF_PART_BORDER (these are back-capable)
    for pb in part_borders:
        pb_vertices = list(pb.vertices())
        if not pb_vertices: continue
        pb_bbox = get_bbox(pb_vertices)
        pb_center = get_bbox_center(pb_vertices)

        # Confirm the part border is within the back sheet border
        if pb_center and is_point_inside_bbox(pb_center, back_sheet_bbox, tolerance=1.1):
            # Find matching cutting line border based ONLY on dimensions
            pb_dims = get_bbox_dimensions_sorted(pb_vertices)
            matched_cl_border = None
            matched_cl_index = -1

            # Search for a cutting line with matching dimensions
            for idx, cl_border in enumerate(cutting_lines_borders):
                 if idx not in used_cutting_lines_indices:
                      cl_vertices = list(cl_border.vertices())
                      if not cl_vertices: continue
                      cl_dims = get_bbox_dimensions_sorted(cl_vertices)
                      
                      # Check only for dimension match, ignore position
                      if math.isclose(pb_dims[0], cl_dims[0], abs_tol=0.1) and \
                         math.isclose(pb_dims[1], cl_dims[1], abs_tol=0.1):
                           matched_cl_border = cl_border
                           matched_cl_index = idx
                           print(f"DEBUG:   تطابق ابعاد پیدا شد برای مرز {part_border_layer} (ابعاد {pb_dims[0]:.1f} x {pb_dims[1]:.1f}) با مرز {cutting_lines_layer} (ابعاد {cl_dims[0]:.1f} x {cl_dims[1]:.1f}).")
                           break

            if matched_cl_border:
                 grouped_panels.append({
                     'borders': [pb, matched_cl_border],
                     'type': 'back_capable',
                     'primary_border': pb,
                     'secondary_border': matched_cl_border,
                     'primary_bbox': get_bbox(list(pb.vertices())),
                     'secondary_bbox': get_bbox(list(matched_cl_border.vertices())),
                     'sheet_border_front_bbox': front_sheet_bbox,
                     'sheet_border_back_bbox': back_sheet_bbox
                 })
                 used_cutting_lines_indices.add(matched_cl_index)
                 print(f"DEBUG: پنل back_capable گروه بندی شد با مرز ABF_CUTTING_LINES هم‌اندازه.")
            else:
                 print(f"WARNING: مرز ABF_PART_BORDER هم‌اندازه‌ای در لایه ABF_CUTTING_LINES یافت نشد.")
                 grouped_panels.append({
                     'borders': [pb],
                     'type': 'back_capable_standalone',
                     'primary_border': pb,
                     'secondary_border': None,
                     'primary_bbox': get_bbox(list(pb.vertices())),
                     'secondary_bbox': None,
                     'sheet_border_front_bbox': front_sheet_bbox,
                     'sheet_border_back_bbox': back_sheet_bbox
                 })
                 print(f"DEBUG: پنل back_capable_standalone گروه بندی شد.")

    # 2. Process remaining _ABF_CUTTING_LINES (front-only panels)
    for idx, cl_border in enumerate(cutting_lines_borders):
        if idx not in used_cutting_lines_indices:
            cl_vertices = list(cl_border.vertices())
            if not cl_vertices: continue
            cl_bbox = get_bbox(cl_vertices)
            cl_center = get_bbox_center(cl_vertices)

            if cl_center and is_point_inside_bbox(cl_center, front_sheet_bbox, tolerance=1.1):
                 grouped_panels.append({
                     'borders': [cl_border],
                     'type': 'front_only',
                     'primary_border': cl_border,
                     'secondary_border': None,
                     'primary_bbox': cl_bbox,
                     'secondary_bbox': None,
                     'sheet_border_front_bbox': front_sheet_bbox,
                     'sheet_border_back_bbox': back_sheet_bbox
                 })
                 print(f"DEBUG: مرز ABF_CUTTING_LINES به عنوان پنل فیزیکی تنها front_only اضافه شد.")
            else:
                 print(f"WARNING: مرز ABF_CUTTING_LINES خارج از محدوده ورق جلو یافت شد.")

    print(f"\nDEBUG: تعداد کل پنل‌های فیزیکی شناسایی شده پس از گروه بندی: {len(grouped_panels)}")
    return grouped_panels


def convert_coords_to_panel_system(entity_x, entity_y, panel_type, panel_xml_length, panel_xml_width, primary_bbox_panel, secondary_bbox_panel, sheet_border_front_bbox, sheet_border_back_bbox, tolerance=0.1):
    """
    Converts DXF coordinates to XML panel coordinates with proper scaling and bounds checking.
    Returns transformed coordinates (rel_x, rel_y) and determined face (5 or 6).
    """
    print(f"DEBUG:   Converting DXF point ({entity_x:.3f}, {entity_y:.3f}) for panel (Type: {panel_type})")
    print(f"DEBUG:     Panel XML dimensions: {panel_xml_length:.0f} x {panel_xml_width:.0f}")

    # Determine which sheet the point is in
    in_front = is_point_inside_bbox((entity_x, entity_y), sheet_border_front_bbox, tolerance)
    in_back = is_point_inside_bbox((entity_x, entity_y), sheet_border_back_bbox, tolerance)
    
    # Select reference bbox and face based on sheet location
    reference_bbox = None
    face = None
    
    if in_front:
        reference_bbox = secondary_bbox_panel if panel_type == 'back_capable' else primary_bbox_panel
        face = "5"
    elif in_back:
        reference_bbox = primary_bbox_panel
        face = "6"
    else:
        print(f"DEBUG:     Point not in any sheet bbox")
        return None, None, None
        
    if not reference_bbox:
        print(f"DEBUG:     No valid reference bbox")
        return None, None, None

    # Get reference bbox dimensions and current point position
    ref_min_x, ref_min_y, ref_max_x, ref_max_y = reference_bbox
    ref_width = ref_max_x - ref_min_x
    ref_height = ref_max_y - ref_min_y
    
    # Calculate relative position within reference bbox (0-1)
    if ref_width == 0 or ref_height == 0:
        print(f"DEBUG:     Invalid reference bbox dimensions")
        return None, None, None
        
    rel_pos_x = (entity_x - ref_min_x) / ref_width
    rel_pos_y = (entity_y - ref_min_y) / ref_height
    
    # Scale to panel dimensions
    xml_x = rel_pos_x * panel_xml_length
    xml_y = rel_pos_y * panel_xml_width
    
    # Mirror Y coordinate for back face
    if face == "6":
        xml_y = panel_xml_width - xml_y
        
    # Ensure coordinates are within bounds
    xml_x = max(0.0, min(xml_x, panel_xml_length))
    xml_y = max(0.0, min(xml_y, panel_xml_width))
    
    print(f"DEBUG:     Transformed to XML coordinates: ({xml_x:.3f}, {xml_y:.3f}), Face: {face}")
    return xml_x, xml_y, face


def create_drilling_xml(machines_element, entity, panel_type, panel_length, panel_width, primary_bbox_panel, secondary_bbox_panel, sheet_border_front_bbox, sheet_border_back_bbox, tolerance, config):
    """
    Processes a CIRCLE entity (ABF_D*) and adds a Machining Type=2 (Drilling) element to XML.
    Determines Face based on which sheet border the entity is within, and Depth from layer name.
    Uses layer names from config.
    Ensures all coordinates are within panel bounds.
    """
    center = entity.dxf.center # Circle center coordinates (includes Z)
    # Use the coordinate conversion function that maps relative to PANEL border bounds and determines Face.
    rel_x, rel_y, face = convert_coords_to_panel_system(center.x, center.y, panel_type, panel_length, panel_width, primary_bbox_panel, secondary_bbox_panel, sheet_border_front_bbox, sheet_border_back_bbox, tolerance)

    if rel_x is None or rel_y is None or face is None:
        print(f"DEBUG:   ایجاد Machining Type=2 (Drilling) برای موجودیت {entity.dxftype()} در لایه {entity.dxf.layer} به دلیل عدم تبدیل مختصات یا Face ناموفق بود.")
        return # Skip if coordinate conversion or face determination failed

    # Use the passed panel_length and panel_width for the debug message
    print(f"DEBUG:   مرکز دایره (DXF): ({center.x:.3f}, {center.y:.3f}, Z={center.z:.3f}) -> تبدیل شده (پنل XML - به سیستم {panel_length:.0f}x{panel_width:.0f}): ({rel_x:.3f}, {rel_y:.3f}), Face: {face}") # Updated Debug message

    layer_name_upper = entity.dxf.layer.upper()
    diameter = round(entity.dxf.radius * 2, 3) # Use the actual circle diameter from DXF

    # Get depth directly from the number after D in the layer name
    depth_number_from_layer = get_number_from_layer_name_after_D(layer_name_upper, config['drilling_prefix'])
    depth = str(depth_number_from_layer) if depth_number_from_layer is not None else "0" # Use the extracted number as depth

    print(f"DEBUG:   پنل {panel_type}، موجودیت {entity.dxftype()} روی لایه {entity.dxf.layer} (قطر {diameter}) -> Face {face}, عمق از لایه: {depth}") # Use extracted depth in debug

    print(f"DEBUG:   در حال ایجاد Machining Type=2 (Drilling) در XML برای X={rel_x:.3f}, Y={rel_y:.3f}, Face={face}, Diameter={diameter:.3f}, Depth={depth}")
    ET.SubElement(machines_element, "Machining",
                  Type="2", # Drilling
                  IsGenCode="2",
                  Face=face, # Use the determined face
                  X=f"{rel_x:.3f}",
                  Y=f"{rel_y:.3f}",
                  Diameter=f"{diameter:.3f}", # Use actual diameter from DXF entity
                  Depth=depth) # Use depth derived from layer name


def create_pocket_xml(machines_element, entity, panel_length, panel_width, panel_type, primary_bbox_panel, secondary_bbox_panel, sheet_border_front_bbox, sheet_border_back_bbox, tolerance, config):
    """
    Processes an ABF_DSIDE_8 LWPOLYLINE entity and adds a Machining Type=1 (Pocket/Routing) element to XML.
    Determines Face based on proximity to panel edges (Faces 1, 2, 3, 4) in the XML system.
    Uses layer names from config.
    """
    vertices = list(entity.vertices())
    # Ensure it's a closed polyline with at least 4 vertices
    if len(vertices) < 4 or not entity.dxf.flags & 1:
        print(f"DEBUG:   LWPOLYLINE روی لایه {config['pocket_layer']} یک مستطیل بسته معتبر نیست یا تعداد نقاط کافی ندارد. رد می شود.")
        return

    rect_min_x = rect_min_y = float('inf')
    rect_max_x = rect_max_y = -float('inf')
    for v in vertices:
        rect_min_x = min(rect_min_x, v[0])
        rect_min_y = min(rect_min_y, v[1])
        rect_max_x = max(rect_max_x, v[0])
        rect_max_y = max(rect_max_y, v[1])

    pocket_z_xml = "8" # Z height from bottom of panel (fixed based on layer name)
    rect_dx = rect_max_x - rect_min_x
    rect_dy = rect_max_y - rect_min_y
    # Pocket depth in XML is the larger dimension of the DXF rectangle
    pocket_depth_xml = max(rect_dx, rect_dy)
    pocket_diameter_xml = "8" # Pocket tool diameter in XML is usually a fixed value (e.g., 8 for ABF_DSIDE_8)

    center_orig_x = (rect_min_x + rect_max_x) / 2.0
    center_orig_y = (rect_min_y + rect_max_y) / 2.0
    # Convert Pocket center to local panel coordinate system (Length x Width)
    # Note: Use the standard tolerance for coordinate conversion
    # Pass PANEL bboxes and sheet bboxes for conversion
    center_rel_x, center_rel_y, determined_face_from_sheet = convert_coords_to_panel_system(center_orig_x, center_orig_y, panel_type, panel_length, panel_width, primary_bbox_panel, secondary_bbox_panel, sheet_border_front_bbox, sheet_border_back_bbox, tolerance)

    if center_rel_x is None or center_rel_y is None or determined_face_from_sheet is None:
        print(f"DEBUG:   ایجاد Machining Type=1 (Pocket) برای موجودیت {entity.dxftype()} در لایه {entity.dxf.layer} به دلیل عدم تبدیل مختصات یا Face ناموفق بود.")
        return # Skip if coordinate conversion or face determination failed

    print(f"DEBUG:   مرکز Pocket (DXF): ({center_orig_x:.3f}, {center_orig_y:.3f}) -> تبدیل شده (پنل XML - به سیستم {panel_length:.0f}x{panel_width:.0f}): ({center_rel_x:.3f}, {center_rel_y:.3f}), Face از ورق: {determined_face_from_sheet}") # Updated Debug message

    # --- Determine Face (1, 2, 3, or 4) for Type=1 based on position in XML panel system ---
    # As per user's strict rule: Type 1 must have Faces 1, 2, 3, or 4.
    # Determine Face based on proximity to XML panel edges using an increased tolerance.
    face = None # Start with no face assigned

    # Use a larger tolerance for determining which edge the pocket is on
    pocket_edge_tolerance = 30.0 # Increased tolerance (adjust if needed based on DXF layout)

    is_close_to_min_y_panel = math.isclose(center_rel_y, 0, abs_tol=pocket_edge_tolerance) # Near Y=0 (Bottom)
    is_close_to_max_y_panel = math.isclose(center_rel_y, panel_width, abs_tol=pocket_edge_tolerance) # Near Y=Width (Top)
    is_close_to_min_x_panel = math.isclose(center_rel_x, 0, abs_tol=pocket_edge_tolerance) # Near X=0 (Left)
    is_close_to_max_x_panel = math.isclose(center_rel_x, panel_length, abs_tol=pocket_edge_tolerance) # Near X=Length (Right)

    print(f"DEBUG:   بررسی نزدیکی Pocket به لبه های پنل در سیستم پنل (تلرانس {pocket_edge_tolerance}): min_y_panel={is_close_to_min_y_panel}, max_y_panel={is_close_to_max_y_panel}, min_x_panel={is_close_to_min_x_panel}, max_x_panel={is_close_to_max_x_panel}")

    # Assign Face based on proximity to XML edges (Bottom->2, Top->1, Right->3, Left->4)
    # Prioritize Y-edges (Bottom/Top) if close to a corner (adjust priority if needed)
    # Based on the sample XML and common practice:
    # Face 1: Top edge (Y=Width in XML)
    # Face 2: Bottom edge (Y=0 in XML)
    # Face 3: Right edge (X=Length in XML)
    # Face 4: Left edge (X=0 in XML)

    if is_close_to_max_y_panel: face = "1" # Top Edge (Y=Width)
    elif is_close_to_min_y_panel: face = "2" # Bottom Edge (Y=0)
    elif is_close_to_max_x_panel: face = "3" # Right Edge (X=Length)
    elif is_close_to_min_x_panel: face = "4" # Left Edge

    print(f"DEBUG:   Assigned Face after proximity check: {face}") # Added debug print


    if face is None:
         # If the Pocket is not close to any edge within the increased tolerance.
         # Assign a default Face (e.g., "5") and log a warning.
         # As per user's requirement, Type=1 MUST have Face 1-4.
         # If no edge is detected, this indicates a potential issue in the DXF or tolerance.
         # Let's assign Face 5 and keep the warning, as forcing 1-4 might be incorrect for an internal pocket.
         # However, if these are truly edge pockets that aren't detected, the tolerance might need increasing further.
         face = "5" # Default Face 5 if not close to any edge
         print(f"DEBUG:   هشدار: Pocket در مختصات تبدیل شده ({center_rel_x:.3f}, {center_rel_y:.3f}) به هیچ لبه پنل نزدیک نیست (با تلرانس {pocket_edge_tolerance}). Face پیش فرض: {face} (5).")


    # X and Y coordinates for Pockets in XML should be on the edge for Faces 1-4.
    # Adjust X or Y based on the determined Face.
    xml_x_val = center_rel_x
    xml_y_val = center_rel_y

    # Snap the XML coordinate to the edge coordinate based on the assigned face (if Face is 1-4)
    if face == "2": # Bottom Edge (Y=0 in XML)
        xml_y_val = 0.0
    elif face == "1": # Top Edge (Y=Width in XML)
        xml_y_val = panel_width
    elif face == "3": # Right Edge (X=Length in XML)
         xml_x_val = panel_length
    elif face == "4": # Left Edge (X=0 in XML)
         xml_x_val = 0.0
    # Note: If face is 5, xml_x_val and xml_y_val remain center_rel_x and center_rel_y.

    print(f"DEBUG:   در حال ایجاد Machining Type=1 (Pocket) در XML برای X={xml_x_val:.3f}, Y={xml_y_val:.3f}, Z={pocket_z_xml}, Face={face}, Diameter={float(pocket_diameter_xml):.3f}, Depth={pocket_depth_xml:.3f}")
    ET.SubElement(machines_element, "Machining",
                      Type="1", # Pocket/Routing
                      IsGenCode="2",
                      Face=face, # Use the determined face (1-4 or 5)
                      X=f"{xml_x_val:.3f}",
                      Y=f"{xml_y_val:.3f}",
                      Z=pocket_z_xml,
                      Diameter=f"{float(pocket_diameter_xml):.3f}", # Tool diameter, format as float
                      Depth=f"{pocket_depth_xml:.3f}") # Cutting depth


def create_groove_xml(machines_element, entity, panel_length, panel_width, panel_type, primary_bbox_panel, secondary_bbox_panel, sheet_border_front_bbox, sheet_border_back_bbox, tolerance, panel_thickness, config):
    """
    Processes an ABF_GROOVE LWPOLYLINE entity and adds a Machining Type=4 (Groove) element to XML.
    Determines the groove centerline from the LWPOLYLINE bounding box.
    Uses layer names from config. Face is typically 5 (front).
    """
    vertices = list(entity.vertices())
    # Ensure it's a closed polyline with at least 4 vertices (assuming rectangular groove representation)
    if len(vertices) < 4 or not entity.dxf.flags & 1:
        print(f"DEBUG:   LWPOLYLINE روی لایه {config['groove_layer']} یک مستطیل بسته معتبر نیست یا تعداد نقاط کافی ندارد. رد می شود.")
        return

    # Calculate the bounding box of the LWPOLYLINE in DXF
    rect_min_x, rect_min_y, rect_max_x, rect_max_y = get_bbox(vertices)

    rect_dx = rect_max_x - rect_min_x
    rect_dy = rect_max_y - rect_min_y

    # Calculate groove width from the shorter dimension of the bounding box
    groove_width = min(rect_dx, rect_dy)
    half_groove_width = groove_width / 2.0

    # Determine the start and end points of the groove centerline in DXF coordinates
    start_point_dxf = None
    end_point_dxf = None

    # Use the corners of the bounding box that define the longer side,
    # but adjust the coordinate perpendicular to the groove direction by half groove width
    if rect_dy > rect_dx: # Longer dimension is along DXF Y (Groove is vertically oriented in DXF)
        # Centerline is along DXF X at rect_min_x + half_groove_width
        # Start and end points are at rect_min_y and rect_max_y
        start_point_dxf = (rect_min_x + half_groove_width, rect_min_y)
        end_point_dxf = (rect_min_x + half_groove_width, rect_max_y)
        print(f"DEBUG:   شیار LWPOLYLINE در DXF عمودی است. نقاط خط مرکزی استفاده شده: ({start_point_dxf[0]:.3f}, {start_point_dxf[1]:.3f}) تا ({end_point_dxf[0]:.3f}, {end_point_dxf[1]:.3f})")
    else: # Longer dimension is along DXF X (or it's a square) (Groove is horizontally oriented in DXF)
        # Centerline is along DXF Y at rect_min_y + half_groove_width
        # Start and end points are at rect_min_x and rect_max_x
        start_point_dxf = (rect_min_x, rect_min_y + half_groove_width)
        end_point_dxf = (rect_max_x, rect_min_y + half_groove_width)
        print(f"DEBUG:   شیار LWPOLYLINE در DXF افقی است. نقاط خط مرکزی استفاده شده: ({start_point_dxf[0]:.3f}, {start_point_dxf[1]:.3f}) تا ({end_point_dxf[0]:.3f}, {end_point_dxf[1]:.3f})")


    if start_point_dxf is None or end_point_dxf is None:
         print(f"DEBUG:   نمی توان نقاط شروع/پایان شیار را از LWPOLYLINE در لایه {config['groove_layer']} تعیین کرد. رد می شود.")
         return # Should not happen for a valid rectangle, but as a safeguard


    # Convert start and end points from DXF (centerline) to XML panel system using the new logic
    # Pass PANEL bboxes and sheet bboxes for conversion
    start_x_xml, start_y_xml, face_start = convert_coords_to_panel_system(start_point_dxf[0], start_point_dxf[1], panel_type, panel_length, panel_width, primary_bbox_panel, secondary_bbox_panel, sheet_border_front_bbox, sheet_border_back_bbox, tolerance=1.0)
    end_x_xml, end_y_xml, face_end = convert_coords_to_panel_system(end_point_dxf[0], end_point_dxf[1], panel_type, panel_length, panel_width, primary_bbox_panel, secondary_bbox_panel, sheet_border_front_bbox, sheet_border_back_bbox, tolerance=1.0)


    if start_x_xml is None or start_y_xml is None or end_x_xml is None or end_y_xml is None or face_start is None or face_end is None:
        print(f"DEBUG:   تبدیل مختصات برای نقاط شروع/پایان شیار در لایه {config['groove_layer']} ناموفق بود. رد می شود.")
        # Add more specific debug here to see which point failed conversion
        # print(f"DEBUG:     Start point DXF: ({start_point_dxf[0]:.3f}, {start_point_dxf[1]:.3f}) -> XML: ({start_x_xml}, {start_y_xml})")
        # print(f"DEBUG:     End point DXF: ({end_point_dxf[0]:.3f}, {end_point_dxf[1]:.3f}) -> XML: ({end_x_xml}, {end_y_xml})")
        return # Skip if coordinate conversion failed

    # For grooves, assume a single face for the entire groove.
    # If start and end points are in different sheet bboxes, this might indicate an issue or a groove spanning faces.
    # For simplicity, let's assume grooves are entirely on Face 5 or Face 6 based on the start point's sheet bbox.
    groove_face = face_start # Use the face determined for the start point

    print(f"DEBUG:   شیار (DXF - خط مرکزی): ({start_point_dxf[0]:.3f}, {start_point_dxf[1]:.3f}) تا ({end_point_dxf[0]:.3f}, {end_point_dxf[1]:.3f}) -> تبدیل شده (پنل XML - به سیستم {panel_length:.0f}x{panel_width:.0f}): ({start_x_xml:.3f}, {start_y_xml:.3f}) تا ({end_x_xml:.3f}, {end_y_xml:.3f})")


    # Fixed attributes based on sample XML and layer type
    groove_type = "4"
    is_gen_code = "2"
    # face is determined by sheet bbox now (groove_face)
    z_height = f"{panel_thickness:.0f}" # Z height often corresponds to material thickness for grooves
    end_z_height = f"{panel_thickness:.0f}" # EndZ is panel thickness
    depth = "8" # Fixed depth from sample (This could potentially come from layer name or other DXF attribute)
    tool_offset = "中" # Fixed from sample

    # Determine Groove width for XML - should be the calculated width
    xml_groove_width = f"{groove_width:.3f}" # Format as float

    print(f"DEBUG:   در حال ایجاد Machining Type=4 (Groove) در XML برای X={start_x_xml:.3f}, Y={start_y_xml:.3f}, EndX={end_x_xml:.3f}, EndY={end_y_xml:.3f}, Face={groove_face}, Width={xml_groove_width}")
    ET.SubElement(machines_element, "Machining",
                  Type=groove_type,
                  IsGenCode=is_gen_code,
                  Face=groove_face, # Use the determined face
                  X=f"{start_x_xml:.3f}",
                  Y=f"{start_y_xml:.3f}",
                  Z=z_height,
                  EndX=f"{end_x_xml:.3f}",
                  EndY=f"{end_y_xml:.3f}",
                  EndZ=end_z_height,
                  Depth=depth,
                  Width=xml_groove_width, # Use calculated groove width
                  ToolOffset=tool_offset)


def process_machining_entities_for_panel(doc, panel_element, panel_group_info, panel_length, panel_width, panel_thickness, config):
    """
    Iterates through all entities in the DXF model space, identifies
    those belonging to the current panel group, and processes them
    into Machining XML elements by calling appropriate helper functions.
    Uses layer names from config.
    """
    machines_element = panel_element.find('Machines') # Get the already created Machines element
    borders_in_group = panel_group_info['borders']
    panel_type = panel_group_info['type']
    # Get sheet border bboxes from panel_group_info
    sheet_border_front_bbox = panel_group_info['sheet_border_front_bbox']
    sheet_border_back_bbox = panel_group_info['sheet_border_back_bbox']
    # Get PANEL bboxes from panel_group_info
    primary_bbox_panel = panel_group_info['primary_bbox']
    secondary_bbox_panel = panel_group_info.get('secondary_bbox')

    # Calculate overall BBox for all borders in this panel group (for entity containment check)
    group_min_x, group_min_y = float('inf'), float('inf')
    group_max_x, group_max_y = -float('inf'), -float('inf')

    for border_entity_in_group in borders_in_group:
         border_vertices = list(border_entity_in_group.vertices())
         if not border_vertices: continue # Skip empty borders
         min_x, min_y, max_x, max_y = get_bbox(border_vertices)
         group_min_x = min(group_min_x, min_x)
         group_min_y = min(group_min_y, min_y)
         group_max_x = max(group_max_x, max_x)
         group_max_y = max(group_max_y, max_y)

    # Use the tolerance from convert_coords_to_panel_system (default 1.0) for containment check
    tolerance = 1.0 # Explicitly use the conversion tolerance here

    print(f"DEBUG: اسکان و پردازش موجودیت‌های ماشینکاری برای پنل فیزیکی...")

    for entity in doc.modelspace():
        # Debug message for every entity read from DXF
        if hasattr(entity, 'dxf') and hasattr(entity.dxf, 'layer'):
            print(f"DEBUG:   در حال بررسی موجودیت DXF: نوع={entity.dxftype()}, لایه={entity.dxf.layer}")
        else:
            print(f"DEBUG:   در حال بررسی موجودیت DXF: نوع={entity.dxftype()} (بدون اطلاعات لایه)")

        # Skip border entities themselves and sheet borders
        if entity in borders_in_group or entity.dxftype() == 'LWPOLYLINE' and entity.dxf.layer.upper() == config['sheet_border'].upper():
            print(f"DEBUG:     موجودیت {entity.dxftype()} در لایه {entity.dxf.layer} یک مرز پنل یا مرز ورق است. رد می شود.")
            continue

        # Get the entity's reference point for checking containment within the panel group's overall BBox
        entity_point_to_check = None
        if hasattr(entity, 'dxf'):
             try:
                 if entity.dxftype() == 'CIRCLE':
                     entity_point_to_check = (entity.dxf.center.x, entity.dxf.center.y)
                 elif entity.dxftype() == 'LWPOLYLINE':
                     vertices = list(entity.vertices())
                     if vertices:
                         # Use the center of the polyline's bounding box
                         min_x, min_y, max_x, max_y = get_bbox(vertices)
                         entity_point_to_check = ((min_x + max_x) / 2.0, (min_y + max_y) / 2.0)
                 elif entity.dxftype() == 'LINE':
                     # Use midpoint of line
                     entity_point_to_check = ((entity.dxf.start.x + entity.dxf.end.x) / 2.0,
                                            (entity.dxf.start.y + entity.dxf.end.y) / 2.0)
                 elif entity.dxftype() in ['ARC', 'ELLIPSE']:
                     # Use center point
                     entity_point_to_check = (entity.dxf.center.x, entity.dxf.center.y)
                 elif entity.dxftype() == 'POINT':
                     entity_point_to_check = (entity.dxf.location.x, entity.dxf.location.y)
                 else:
                     print(f"DEBUG:     نوع موجودیت {entity.dxftype()} پشتیبانی نمی شود یا نقطه مرجع برای آن تعریف نشده است.")
             except Exception as e:
                 print(f"DEBUG:     خطا در استخراج نقطه مرجع برای موجودیت {entity.dxftype()}: {str(e)}")

        is_within_panel_group = False
        if entity_point_to_check is not None and len(entity_point_to_check) >= 2:
             # Check if the reference point is within either sheet border bbox
             is_within_front = is_point_inside_bbox(entity_point_to_check, sheet_border_front_bbox, tolerance)
             is_within_back = is_point_inside_bbox(entity_point_to_check, sheet_border_back_bbox, tolerance)
             is_within_panel_group = is_within_front or is_within_back
             print(f"DEBUG:     نقطه مرجع موجودیت: ({entity_point_to_check[0]:.3f}, {entity_point_to_check[1]:.3f})")
             print(f"DEBUG:     درون محدوده ورق جلو: {is_within_front}, درون محدوده ورق پشت: {is_within_back}")

        if not is_within_panel_group:
            print(f"DEBUG:     موجودیت {entity.dxftype()} در لایه {entity.dxf.layer} خارج از محدوده پنل است. رد می شود.")
            continue

        # --- Process machining entities (only for entities within the overall BBox) ---

        # Process CIRCLEs (any layer starting with ABF_D) - Type="2".
        if entity.dxftype() == 'CIRCLE' and entity.dxf.layer.upper().startswith(config['drilling_prefix'].upper()):
            create_drilling_xml(machines_element, entity, panel_type, panel_length, panel_width, 
                              primary_bbox_panel, secondary_bbox_panel, sheet_border_front_bbox, 
                              sheet_border_back_bbox, tolerance, config)
        elif entity.dxftype() == 'LWPOLYLINE' and entity.dxf.layer.upper() == config['pocket_layer'].upper():
            create_pocket_xml(machines_element, entity, panel_length, panel_width, panel_type,
                            primary_bbox_panel, secondary_bbox_panel, sheet_border_front_bbox,
                            sheet_border_back_bbox, tolerance, config)
        elif entity.dxftype() == 'LWPOLYLINE' and entity.dxf.layer.upper() == config['groove_layer'].upper():
            create_groove_xml(machines_element, entity, panel_length, panel_width, panel_type,
                            primary_bbox_panel, secondary_bbox_panel, sheet_border_front_bbox,
                            sheet_border_back_bbox, tolerance, panel_thickness, config)
        else:
            print(f"DEBUG:     موجودیت {entity.dxftype()} در لایه {entity.dxf.layer} با الگوهای ماشینکاری شناخته شده مطابقت ندارد. رد می شود.")


def create_panel_xml_structure(panel_id, panel_name, length, width, thickness):
    """
    Creates the basic XML structure for a panel with given dimensions.
    Returns the root 'Panel' element and its 'Machines' child element.
    """
    # Create root Panel element with attributes
    panel_element = ET.Element("Panel",
                             ID=str(panel_id),
                             Name=panel_name,
                             Length=f"{length:.3f}",
                             Width=f"{width:.3f}",
                             Thickness=f"{thickness:.3f}",
                             CreationDate="",
                             Grain="0",
                             MaterialType="0")

    # Add required child elements
    # Add Machines element (will contain machining operations)
    machines_element = ET.SubElement(panel_element, "Machines")
    
    # Add EmptySlots element (required but typically empty)
    empty_slots = ET.SubElement(panel_element, "EmptySlots")
    
    # Add Parameters element (required but typically empty)
    parameters = ET.SubElement(panel_element, "Parameters")

    return panel_element, machines_element


# --- Main Conversion Function ---

def dxf_to_custom_xml(input_file, config, panel_thickness=16.0): # Removed unused output_base_name
    """
    Main function to read DXF file, identify and process panels and their
    machining entities, and generate corresponding XML files.
    Uses layer names from config.
    """
    try:
        # Load the DXF document
        doc = ezdxf.readfile(input_file)
        print(f"DEBUG: فایل DXF '{input_file}' با موفقیت بارگذاری شد.")

        # 1. Find and group physical panels based on sheet borders and part borders
        grouped_panels = find_and_group_panels(doc, config)

        if not grouped_panels:
             print(f"❌ خطا: هیچ پنل فیزیکی برای پردازش یافت نشد.")
             return False # Indicate failure

        # Get the base name of the input DXF file without extension
        dxf_base_name = os.path.splitext(os.path.basename(input_file))[0]
        
        # Create output folder based on DXF file name
        output_folder = os.path.join(os.path.dirname(input_file), dxf_base_name)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
            print(f"DEBUG: Created output folder: {output_folder}")

        # 2. Process each grouped physical panel
        for i, panel_group_info in enumerate(grouped_panels):
            primary_border_entity = panel_group_info['primary_border']
            # Determine XML panel dimensions (Length=larger, Width=smaller) from primary border
            # Note: get_bbox_dimensions_sorted returns (smallest, largest)
            # We need to assign the larger dimension to XML Length and smaller to XML Width
            panel_xml_width, panel_xml_length = get_bbox_dimensions_sorted(list(primary_border_entity.vertices()))
            length = panel_xml_length # Assign larger dimension to Length for XML
            width = panel_xml_width   # Assign smaller dimension to Width for XML

            print(f"\nDEBUG: --- شروع پردازش پنل فیزیکی شماره {i+1} (نوع: {panel_group_info['type']}) ---")
            print(f"DEBUG:   ابعاد پنل در XML (Length x Width): {length:.3f} x {width:.3f}")

            # 3. Create basic XML structure for the panel with the new naming convention
            # File name should be DXF_NAME.LengthxWidth.PanelIndex.xml for uniqueness
            output_file_name_base = f"{dxf_base_name}.{length:.0f}x{width:.0f}.{i+1}"
            output_file = os.path.join(output_folder, f"{output_file_name_base}.xml")

            # Panel ID and Name should be the same as the file name base (without .xml)
            panel_id = output_file_name_base
            panel_name = output_file_name_base


            root, panel_element = create_panel_xml_structure(panel_id, panel_name, length, width, panel_thickness)
            print(f"DEBUG:   ساختار اصلی XML برای پنل '{panel_name}' ایجاد شد.")

            # 4. Process machining entities within this panel group and add them to the XML structure
            # Pass all necessary info including panel bboxes and sheet bboxes
            process_machining_entities_for_panel(doc, panel_element, panel_group_info, length, width, panel_thickness, config)

            # 5. Save the generated XML file
            tree = ET.ElementTree(root)
            # Use indent for pretty printing the XML
            # Use minidom for pretty printing for compatibility with older Python versions
            xml_string = ET.tostring(root, encoding='utf-8').decode('utf-8')
            dom = xml.dom.minidom.parseString(xml_string)
            # Ensure pretty_xml_string is bytes before writing to file in binary mode
            pretty_xml_string = dom.toprettyxml(indent="  ", encoding='utf-8')

            with open(output_file, "wb") as f:
                f.write(pretty_xml_string)


            print(f"✅ فایل '{output_file}' با موفقیت برای پنل فیزیکی شماره {i+1} (نوع: {panel_group_info['type']}) ایجاد شد.")
            print(f"DEBUG: --- پایان پردازش پنل فیزیکی شماره {i+1} ---")

        return True # Indicate success


    except FileNotFoundError:
        print(f"❌ خطا: فایل ورودی '{input_file}' یافت نشد.")
        return False # Indicate failure
    except ezdxf.DXFStructureError:
         print(f"❌ خطا: فایل '{input_file}' یک فایل DXF معتبر نیست یا خراب است.")
         return False # Indicate failure
    except Exception as e:
        print(f"❌ خطا در پردازش فایل DXF: {str(e)}")
        # Print full traceback for easier debugging of unexpected errors
        import traceback
        traceback.print_exc()
        return False # Indicate failure


# --- Terminal UI Class ---
class TerminalUI:
    """Handles user interaction in the terminal."""
    def __init__(self, config):
        self.config = config

    def run(self):
        """Runs the terminal UI for selecting and converting a single DXF file."""
        while True: # Keep the UI running until explicitly exited
            try:
                subprocess.run(['clear' if os.name != 'nt' else 'cls'], shell=True)
            except:
                print('\n' * 100)
                
            print("=========================================")
            print("  مبدل DXF به XML برای پنل‌های چوبی")
            print("=========================================")

            dxf_files = self._get_dxf_files_in_current_directory()

            if not dxf_files:
                print("\n❌  هیچ فایل DXF در پوشه جاری یافت نشد.")
                print("لطفاً فایل‌های DXF را در کنار اسکریپت قرار دهید.")
                input("\nبرای خروج، Enter را فشار دهید...")
                sys.exit()

            print("\nفایل‌های DXF موجود در پوشه جاری:")
            for i, file_name in enumerate(dxf_files):
                print(f"{i + 1}. {file_name}")

            print("\n0. خروج")

            try:
                choice = int(input("\nلطفاً شماره فایل مورد نظر را وارد کنید (0-" + str(len(dxf_files)) + "): "))
                if 0 <= choice <= len(dxf_files):
                    if choice == 0:
                        print("\nخروج از برنامه.")
                        sys.exit()
                    else:
                        selected_file = dxf_files[choice - 1]
                        print(f"\n✅  فایل انتخابی: {selected_file}")
                        print("\nدر حال پردازش...")

                        # Run the conversion
                        success = dxf_to_custom_xml(selected_file, self.config)

                        if success:
                            print("\n✅  فرآیند تبدیل با موفقیت به پایان رسید.")
                        else:
                            print("\n❌  فرآیند تبدیل با مشکل مواجه شد.")

                        input("\nبرای ادامه، Enter را فشار دهید...")
                else:
                    print("⚠️  شماره نامعتبر.")
                    input("\nبرای ادامه، Enter را فشار دهید...")
            except ValueError:
                print("⚠️  ورودی نامعتبر.")
                input("\nبرای ادامه، Enter را فشار دهید...")

    def _get_dxf_files_in_current_directory(self):
        """Lists all files with .dxf extension in the current directory."""
        return [f for f in os.listdir('.') if os.path.isfile(f) and f.lower().endswith('.dxf')]
