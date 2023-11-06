def in_rb(
        id,
        rb_regions
):
    in_rb = False
    for i in range(len(rb_regions)):
        start = rb_regions[i][0]
        end = rb_regions[i][1]

        if (id >= start) and (id < end):
            in_rb = True

    return in_rb


def get_flex_regions(
        rb_regions,
        n_res
):
    flex_regions = list()
    i = 0
    while i < n_res:
        if not in_rb(i, rb_regions):
            start = i
            j = i+1
            while j < n_res and not in_rb(j, rb_regions):
                j = j+1
            end = j
            flex_regions.append([start, end])
            i = j+1
        i = i+1

    return flex_regions