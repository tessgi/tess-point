from tesspoint import TESSPoint, footprint


def test_tesspoint():
    tp = TESSPoint(1, 1, 1)
    tp.pix2radec(footprint())
