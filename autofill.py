from embroider import param, EmbroideryElement
import PyEmb

class RegionSplitter(object):
    """Split a fill region into FillSection objects."""

    class _Row(object):
        """Describes the intersection of the fill region with a line."""

        def __init__(self, segments, row_num):
            """Initialize a _Row

            Parameters:

            segments - an iterable of (Point, Point) tuples describing the
                       start and end of each segment in the intersection.  The
                       Points must all be colinear.
            row_num - an integer row index, used to group segments into FillSections.
            """

            self.sections = sections
            self.row_num = row_num

    def __init__(self, shape, angle, row_spacing):
        """Initialize a RegionSplitter.

        Parameters:

        shape -- a shapely MultiPolygon describing the fill region
        angle -- angle of rows of stitches
        row_spacing -- space between each row of stitches
        """

        self.shape = shape
        self.angle = angle
        self.row_spacing = row_spacing

    def split():
        """Return an iterable of FillSection objects covering the region."""
        pass

class FillSection(object):
    """An subsection of a fill region that can be stitched all in one go."""

    def __init__(self, rows, row_num):
       """Initialize a FillSection

       Parameters:

       rows -- an iterable of (PyEmb.Point, PyEmb.Point) tuples describing the
               start and endpoint of the rows in this subsection
       row_num -- an integer index of the first row, used to determine whether
                  FillSections are neighbors
       """

       self.rows = rows
       self.row_num = row_num

class AutoFillV2(EmbroideryElement):
    def __init__(self, *args, **kwargs):
        super(Fill, self).__init__(*args, **kwargs)

    @property
    @param('auto_fill', 'Manually routed fill stitching', type='toggle', inverse=True, default=True)
    def auto_fill(self):
        return self.get_boolean_param('auto_fill', True)

    @property
    @param('angle', 'Angle of lines of stitches', unit='deg', type='float')
    @cache
    def angle(self):
        return math.radians(self.get_float_param('angle', 0))

    @property
    def color(self):
        return self.get_style("fill")

    @property
    @param('flip', 'Flip fill (start right-to-left)', type='boolean')
    def flip(self):
        return self.get_boolean_param("flip", False)

    @property
    @param('row_spacing_mm', 'Spacing between rows', unit='mm', type='float')
    def row_spacing(self):
        return self.get_float_param("row_spacing_mm")

    @property
    @param('max_stitch_length_mm', 'Maximum fill stitch length', unit='mm', type='float')
    def max_stitch_length(self):
        return self.get_float_param("max_stitch_length_mm")

    @property
    @param('staggers', 'Stagger rows this many times before repeating', type='int')
    def staggers(self):
        return self.get_int_param("staggers", 4)

    @property
    @cache
    def paths(self):
        return self.flatten(self.parse_path())

    @property
    @cache
    def shape(self):
        poly_ary = []
        for sub_path in self.paths:
            point_ary = []
            last_pt = None
            for pt in sub_path:
                if (last_pt is not None):
                    vp = (pt[0] - last_pt[0], pt[1] - last_pt[1])
                    dp = math.sqrt(math.pow(vp[0], 2.0) + math.pow(vp[1], 2.0))
                    # dbg.write("dp %s\n" % dp)
                    if (dp > 0.01):
                        # I think too-close points confuse shapely.
                        point_ary.append(pt)
                        last_pt = pt
                else:
                    last_pt = pt
            poly_ary.append(point_ary)

        # shapely's idea of "holes" are to subtract everything in the second set
        # from the first. So let's at least make sure the "first" thing is the
        # biggest path.
        # TODO: actually figure out which things are holes and which are shells
        poly_ary.sort(key=lambda point_list: shgeo.Polygon(point_list).area, reverse=True)

        polygon = shgeo.MultiPolygon([(poly_ary[0], poly_ary[1:])])
        # print >> sys.stderr, "polygon valid:", polygon.is_valid
        return polygon

    @cache
    def east(self, angle):
        # "east" is the name of the direction that is to the right along a row
        return PyEmb.Point(1, 0).rotate(-angle)

    @cache
    def north(self, angle):
        return self.east(angle).rotate(math.pi / 2)

    def adjust_stagger(self, stitch, angle, row_spacing, max_stitch_length):
        row_num = round((stitch * self.north(angle)) / row_spacing)
        row_stagger = row_num % self.staggers
        stagger_offset = (float(row_stagger) / self.staggers) * max_stitch_length
        offset = ((stitch * self.east(angle)) - stagger_offset) % max_stitch_length

        return stitch - offset * self.east(angle)

    def intersect_region_with_grating(self, angle=None, row_spacing=None):
        if angle is None:
            angle = self.angle

        if row_spacing is None:
            row_spacing = self.row_spacing

        # the max line length I'll need to intersect the whole shape is the diagonal
        (minx, miny, maxx, maxy) = self.shape.bounds
        upper_left = PyEmb.Point(minx, miny)
        lower_right = PyEmb.Point(maxx, maxy)
        length = (upper_left - lower_right).length()
        half_length = length / 2.0

        # Now get a unit vector rotated to the requested angle.  I use -angle
        # because shapely rotates clockwise, but my geometry textbooks taught
        # me to consider angles as counter-clockwise from the X axis.
        direction = PyEmb.Point(1, 0).rotate(-angle)

        # and get a normal vector
        normal = direction.rotate(math.pi / 2)

        # I'll start from the center, move in the normal direction some amount,
        # and then walk left and right half_length in each direction to create
        # a line segment in the grating.
        center = PyEmb.Point((minx + maxx) / 2.0, (miny + maxy) / 2.0)

        # I need to figure out how far I need to go along the normal to get to
        # the edge of the shape.  To do that, I'll rotate the bounding box
        # angle degrees clockwise and ask for the new bounding box.  The max
        # and min y tell me how far to go.

        _, start, _, end = affinity.rotate(self.shape, angle, origin='center', use_radians=True).bounds

        # convert start and end to be relative to center (simplifies things later)
        start -= center.y
        end -= center.y

        # offset start slightly so that rows are always an even multiple of
        # row_spacing_px from the origin.  This makes it so that abutting
        # fill regions at the same angle and spacing always line up nicely.
        start -= (start + normal * center) % row_spacing

        rows = []

        while start < end:
            p0 = center + normal * start + direction * half_length
            p1 = center + normal * start - direction * half_length
            grating_line = shgeo.LineString((p0, p1))

            res = grating_line.intersection(self.shape)

            if (isinstance(res, shgeo.MultiLineString)):
                runs = map(lambda line_string: line_string.coords, res.geoms)
            else:
                if res.is_empty or len(res.coords) == 1:
                    # ignore if we intersected at a single point or no points
                    start += row_spacing
                    continue
                runs = [res.coords]

            runs = [tuple(PyEmb.Point(x, y) for x, y in run) for run in runs]

            runs.sort(key=lambda seg: seg[0].distance(upper_left))

            if self.flip:
                runs.reverse()
                runs = map(lambda run: tuple(reversed(run)), runs)

            rows.append(runs)

            start += row_spacing

        return rows

    def make_quadrilateral(self, segment1, segment2):
        return shgeo.Polygon((segment1[0], segment1[1], segment2[1], segment2[0], segment1[0]))

    def is_same_run(self, segment1, segment2):
        if shgeo.LineString(segment1).distance(shgeo.LineString(segment2)) > self.row_spacing * 1.1:
            return False

        quad = self.make_quadrilateral(segment1, segment2)
        quad_area = quad.area
        intersection_area = self.shape.intersection(quad).area

        return (intersection_area / quad_area) >= 0.9

    def pull_runs(self, rows):
        # Given a list of rows, each containing a set of line segments,
        # break the area up into contiguous patches of line segments.
        #
        # This is done by repeatedly pulling off the first line segment in
        # each row and calling that a shape.  We have to be careful to make
        # sure that the line segments are part of the same shape.  Consider
        # the letter "H", with an embroidery angle of 45 degrees.  When
        # we get to the bottom of the lower left leg, the next row will jump
        # over to midway up the lower right leg.  We want to stop there and
        # start a new patch.

        # for row in rows:
        #    print >> sys.stderr, len(row)

        # print >>sys.stderr, "\n".join(str(len(row)) for row in rows)

        runs = []
        count = 0
        while (len(rows) > 0):
            run = []
            prev = None

            for row_num in xrange(len(rows)):
                row = rows[row_num]
                first, rest = row[0], row[1:]

                # TODO: only accept actually adjacent rows here
                if prev is not None and not self.is_same_run(prev, first):
                    break

                run.append(first)
                prev = first

                rows[row_num] = rest

            # print >> sys.stderr, len(run)
            runs.append(run)
            rows = [row for row in rows if len(row) > 0]

            count += 1

        return runs

    def section_to_patch(self, group_of_segments, angle=None, row_spacing=None, max_stitch_length=None):
        if max_stitch_length is None:
            max_stitch_length = self.max_stitch_length

        if row_spacing is None:
            row_spacing = self.row_spacing

        if angle is None:
            angle = self.angle

        # print >> sys.stderr, len(groups_of_segments)

        patch = Patch(color=self.color)
        first_segment = True
        swap = False
        last_end = None

        for segment in group_of_segments:
            # We want our stitches to look like this:
            #
            # ---*-----------*-----------
            # ------*-----------*--------
            # ---------*-----------*-----
            # ------------*-----------*--
            # ---*-----------*-----------
            #
            # Each successive row of stitches will be staggered, with
            # num_staggers rows before the pattern repeats.  A value of
            # 4 gives a nice fill while hiding the needle holes.  The
            # first row is offset 0%, the second 25%, the third 50%, and
            # the fourth 75%.
            #
            # Actually, instead of just starting at an offset of 0, we
            # can calculate a row's offset relative to the origin.  This
            # way if we have two abutting fill regions, they'll perfectly
            # tile with each other.  That's important because we often get
            # abutting fill regions from pull_runs().

            (beg, end) = segment

            if (swap):
                (beg, end) = (end, beg)

            row_direction = (end - beg).unit()
            segment_length = beg.distance(end)

            # only stitch the first point if it's a reasonable distance away from the
            # last stitch
            if last_end is None or beg.distance(last_end) > 0.5 * self.options.pixels_per_mm:
                patch.add_stitch(beg)

            first_stitch = self.adjust_stagger(beg, angle, row_spacing, max_stitch_length)

            # we might have chosen our first stitch just outside this row, so move back in
            if (first_stitch - beg) * row_direction < 0:
                first_stitch += row_direction * max_stitch_length

            offset = first_stitch.distance(beg)

            while offset < segment_length:
                patch.add_stitch(beg + offset * row_direction)
                offset += max_stitch_length

            if end.distance(patch.stitches[-1]) > 0.1 * self.options.pixels_per_mm:
                patch.add_stitch(end)

            last_end = end
            swap = not swap

        return patch

    def to_patches(self, last_patch):
        rows_of_segments = self.intersect_region_with_grating()
        groups_of_segments = self.pull_runs(rows_of_segments)

        return [self.section_to_patch(group) for group in groups_of_segments]


class AutoFill(Fill):
    @property
    @param('auto_fill', 'Automatically routed fill stitching', type='toggle', default=True)
    def auto_fill(self):
        return self.get_boolean_param('auto_fill', True)

    @property
    @cache
    def outline(self):
        return self.shape.boundary[0]

    @property
    @cache
    def outline_length(self):
        return self.outline.length

    @property
    def flip(self):
        return False

    @property
    @param('running_stitch_length_mm', 'Running stitch length (traversal between sections)', unit='mm', type='float')
    def running_stitch_length(self):
        return self.get_float_param("running_stitch_length_mm")

    @property
    @param('fill_underlay', 'Underlay', type='toggle', group='AutoFill Underlay')
    def fill_underlay(self):
        return self.get_boolean_param("fill_underlay")

    @property
    @param('fill_underlay_angle', 'Fill angle (default: fill angle + 90 deg)', unit='deg', group='AutoFill Underlay', type='float')
    @cache
    def fill_underlay_angle(self):
        underlay_angle = self.get_float_param("fill_underlay_angle")

        if underlay_angle:
            return math.radians(underlay_angle)
        else:
            return self.angle + math.pi / 2.0

    @property
    @param('fill_underlay_row_spacing_mm', 'Row spacing (default: 3x fill row spacing)', unit='mm', group='AutoFill Underlay', type='float')
    @cache
    def fill_underlay_row_spacing(self):
        return self.get_float_param("fill_underlay_row_spacing_mm") or self.row_spacing * 3

    @property
    @param('fill_underlay_max_stitch_length_mm', 'Max stitch length', unit='mm', group='AutoFill Underlay', type='float')
    @cache
    def fill_underlay_max_stitch_length(self):
        return self.get_float_param("fill_underlay_max_stitch_length_mm" or self.max_stitch_length)

    def validate(self):
        if len(self.shape.boundary) > 1:
            self.fatal("auto-fill: object %s cannot be auto-filled because it has one or more holes.  Please disable auto-fill for this object or break it into separate objects without holes." % self.node.get('id'))

    def is_same_run(self, segment1, segment2):
        if shgeo.Point(segment1[0]).distance(shgeo.Point(segment2[0])) > self.max_stitch_length:
            return False

        if shgeo.Point(segment1[1]).distance(shgeo.Point(segment2[1])) > self.max_stitch_length:
            return False

        return True

    def perimeter_distance(self, p1, p2):
        # how far around the perimeter (and in what direction) do I need to go
        # to get from p1 to p2?

        p1_projection = self.outline.project(shgeo.Point(p1))
        p2_projection = self.outline.project(shgeo.Point(p2))

        distance = p2_projection - p1_projection

        if abs(distance) > self.outline_length / 2.0:
            # if we'd have to go more than halfway around, it's faster to go
            # the other way
            if distance < 0:
                return distance + self.outline_length
            elif distance > 0:
                return distance - self.outline_length
            else:
                # this ought not happen, but just for completeness, return 0 if
                # p1 and p0 are the same point
                return 0
        else:
            return distance

    def connect_points(self, p1, p2):
        patch = Patch(color=self.color)

        pos = self.outline.project(shgeo.Point(p1))
        distance = self.perimeter_distance(p1, p2)
        stitches = abs(int(distance / self.running_stitch_length))

        direction = math.copysign(1.0, distance)
        one_stitch = self.running_stitch_length * direction

        for i in xrange(stitches):
            pos = (pos + one_stitch) % self.outline_length

            stitch = PyEmb.Point(*self.outline.interpolate(pos).coords[0])

            # if we're moving along the fill direction, adjust the stitch to
            # match the fill so it blends in
            if patch.stitches:
                if abs((stitch - patch.stitches[-1]) * self.north(self.angle)) < 0.01:
                    new_stitch = self.adjust_stagger(stitch, self.angle, self.row_spacing, self.max_stitch_length)

                    # don't push the point past the end of this section of the outline
                    if self.outline.distance(shgeo.Point(new_stitch)) <= 0.01:
                        stitch = new_stitch

            patch.add_stitch(stitch)

        return patch

    def get_corner_points(self, section):
        return section[0][0], section[0][-1], section[-1][0], section[-1][-1]

    def nearest_corner(self, section, point):
        return min(self.get_corner_points(section),
                   key=lambda corner: abs(self.perimeter_distance(point, corner)))

    def find_nearest_section(self, sections, point):
        sections_with_nearest_corner = [(i, self.nearest_corner(section, point))
                                        for i, section in enumerate(sections)]
        return min(sections_with_nearest_corner,
                   key=lambda(section, corner): abs(self.perimeter_distance(point, corner)))

    def section_from_corner(self, section, start_corner, angle, row_spacing, max_stitch_length):
        if start_corner not in section[0]:
            section = list(reversed(section))

        if section[0][0] != start_corner:
            section = [list(reversed(row)) for row in section]

        return self.section_to_patch(section, angle, row_spacing, max_stitch_length)

    def do_auto_fill(self, angle, row_spacing, max_stitch_length, starting_point=None):
        rows_of_segments = self.intersect_region_with_grating(angle, row_spacing)
        sections = self.pull_runs(rows_of_segments)

        patches = []
        last_stitch = starting_point
        while sections:
            if last_stitch:
                section_index, start_corner = self.find_nearest_section(sections, last_stitch)
                patches.append(self.connect_points(last_stitch, start_corner))
                patches.append(self.section_from_corner(sections.pop(section_index), start_corner, angle, row_spacing, max_stitch_length))
            else:
                patches.append(self.section_to_patch(sections.pop(0), angle, row_spacing, max_stitch_length))

            last_stitch = patches[-1].stitches[-1]

        return patches

    def to_patches(self, last_patch):
        print >> dbg, "autofill"
        self.validate()

        patches = []

        if last_patch is None:
            last_stitch = None
        else:
            last_stitch = last_patch.stitches[-1]

        if self.fill_underlay:
            patches.extend(self.do_auto_fill(self.fill_underlay_angle, self.fill_underlay_row_spacing, self.fill_underlay_max_stitch_length, last_stitch))
            last_stitch = patches[-1].stitches[-1]

        patches.extend(self.do_auto_fill(self.angle, self.row_spacing, self.max_stitch_length, last_stitch))

        return patches
