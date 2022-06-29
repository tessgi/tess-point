"""Python vectorized version of tess-point"""
import numpy as np
from dataclasses import dataclass
from . import PACKAGEDIR

__all__ = ["TESSPoint", "footprint"]

pointings = {
    key: col
    for col, key in zip(
        np.loadtxt(f"{PACKAGEDIR}/data/pointings.csv", delimiter=",", skiprows=1).T,
        np.loadtxt(
            f"{PACKAGEDIR}/data/pointings.csv", delimiter=",", max_rows=1, dtype=str
        ),
    )
}

tess_params = {
    1: {
        "ang1": 0.101588,
        "ang2": -36.022035,
        "ang3": 90.048315,
        "fl": 145.948116,
        "opt_coef1": 1.00000140,
        "opt_coef2": 0.24779006,
        "opt_coef3": -0.22681254,
        "opt_coef4": 10.78243356,
        "opt_coef5": -34.97817276,
        "x0_ccd1": 31.573417,
        "y0_ccd1": 31.551637,
        "ang_ccd1": 179.980833,
        "x0_ccd2": -0.906060,
        "y0_ccd2": 31.536148,
        "ang_ccd2": 180.000000,
        "x0_ccd3": -31.652818,
        "y0_ccd3": -31.438350,
        "ang_ccd3": -0.024851,
        "x0_ccd4": 0.833161,
        "y0_ccd4": -31.458180,
        "ang_ccd4": 0.001488,
    },
    2: {
        "ang1": -0.179412,
        "ang2": -12.017260,
        "ang3": 90.046500,
        "fl": 145.989933,
        "opt_coef1": 1.00000140,
        "opt_coef2": 0.24069345,
        "opt_coef3": 0.15391120,
        "opt_coef4": 4.05433503,
        "opt_coef5": 3.43136895,
        "x0_ccd1": 31.653635,
        "y0_ccd1": 31.470291,
        "ang_ccd1": 180.010890,
        "x0_ccd2": -0.827405,
        "y0_ccd2": 31.491388,
        "ang_ccd2": 180.000000,
        "x0_ccd3": -31.543794,
        "y0_ccd3": -31.550699,
        "ang_ccd3": -0.006624,
        "x0_ccd4": 0.922834,
        "y0_ccd4": -31.557268,
        "ang_ccd4": -0.015464,
    },
    3: {
        "ang1": 0.066596,
        "ang2": 12.007750,
        "ang3": -89.889085,
        "fl": 146.006602,
        "opt_coef1": 1.00000140,
        "opt_coef2": 0.23452229,
        "opt_coef3": 0.33552009,
        "opt_coef4": 1.92009863,
        "opt_coef5": 12.48880182,
        "x0_ccd1": 31.615486,
        "y0_ccd1": 31.413644,
        "ang_ccd1": 179.993948,
        "x0_ccd2": -0.832993,
        "y0_ccd2": 31.426621,
        "ang_ccd2": 180.000000,
        "x0_ccd3": -31.548296,
        "y0_ccd3": -31.606976,
        "ang_ccd3": 0.000298,
        "x0_ccd4": 0.896018,
        "y0_ccd4": -31.569542,
        "ang_ccd4": -0.006464,
    },
    4: {
        "ang1": 0.030756,
        "ang2": 35.978116,
        "ang3": -89.976802,
        "fl": 146.039793,
        "opt_coef1": 1.00000140,
        "opt_coef2": 0.23920416,
        "opt_coef3": 0.13349450,
        "opt_coef4": 4.77768896,
        "opt_coef5": -1.75114744,
        "x0_ccd1": 31.575820,
        "y0_ccd1": 31.316510,
        "ang_ccd1": 179.968217,
        "x0_ccd2": -0.890877,
        "y0_ccd2": 31.363511,
        "ang_ccd2": 180.000000,
        "x0_ccd3": -31.630470,
        "y0_ccd3": -31.716942,
        "ang_ccd3": -0.024359,
        "x0_ccd4": 0.824159,
        "y0_ccd4": -31.728751,
        "ang_ccd4": -0.024280,
    },
}


@dataclass
class TESSPoint:
    sector: int
    camera: int
    ccd: int

    def __post_init__(self):
        xeul = np.hstack(
            [
                (np.pi / 180.0) * pointings["ra"][pointings["sector"] == self.sector],
                np.pi / 2.0
                - (np.pi / 180.0)
                * pointings["dec"][pointings["sector"] == self.sector],
                (np.pi / 180.0) * pointings["roll"][pointings["sector"] == self.sector]
                + np.pi,
            ]
        )

        self.rmat1 = eulerm323(xeul)
        eulcam = np.asarray(
            [tess_params[self.camera][f"ang{idx}"] for idx in np.arange(1, 4)]
        )
        self.rmat2 = eulerm323(eulcam * (np.pi / 180.0))
        self.rmat4 = np.matmul(self.rmat2, self.rmat1)

    @property
    def opt_coeffs(self):
        return np.asarray(
            [
                tess_params[self.camera][key]
                for key in np.hstack(
                    [["fl"], [f"opt_coef{idx}" for idx in np.arange(1, 6)]]
                )
            ]
        )

    def pix_to_mm(self, coords):
        """convert pixel to mm focal plane position"""
        pixsz = 0.015000
        angle = tess_params[self.camera][f"ang_ccd{self.ccd}"]
        xyb = xyrotate(angle, (coords + 0.5) * pixsz)
        return np.vstack(
            [
                xyb[:, 0] + tess_params[self.camera][f"x0_ccd{self.ccd}"],
                xyb[:, 1] + tess_params[self.camera][f"y0_ccd{self.ccd}"],
            ]
        ).T

    def pix2radec(self, coords):
        xyfp = self.pix_to_mm(coords)
        lng_deg, lat_deg = fp_optics(xyfp, self.opt_coeffs)
        vcam = np.asarray(sphereToCart(lng_deg, lat_deg)).T
        curVec = np.matmul(self.rmat4.T, vcam.T).T
        ra, dec = cartToSphere(curVec)
        return ra / (np.pi / 180.0), dec / (np.pi / 180.0)


def footprint(npoints=50):
    """Gets the column and row points for CCD edges"""
    column = np.hstack(
        [
            np.zeros(npoints),
            np.linspace(0, 2048, npoints),
            np.linspace(0, 2048, npoints),
            np.ones(npoints) * 2048,
        ]
    )
    row = np.hstack(
        [
            np.linspace(0, 2048, npoints),
            np.zeros(npoints),
            np.ones(npoints) * 2048,
            np.linspace(0, 2048, npoints),
        ]
    )
    return np.vstack([column, row]).T


def xyrotate(angle, coords):
    ca = np.cos((np.pi / 180.0) * angle)
    sa = np.sin((np.pi / 180.0) * angle)
    return np.vstack(
        [ca * coords[:, 0] + sa * coords[:, 1], -sa * coords[:, 0] + ca * coords[:, 1]]
    ).T.astype(coords.dtype)


def rev_az_asym(coords):
    asymang = 0.0
    asymfac = 1.0
    xypa = xyrotate(asymang, coords) * np.asarray([1 / asymfac, 1])
    return xyrotate(-asymang, xypa)


def r_of_tanth(z, opt_coeffs):
    tanth = np.tan(z)
    rfp0 = tanth * opt_coeffs[0]
    rfp = np.sum(opt_coeffs[1:] * (tanth ** (2 * np.arange(5))[:, None]).T, axis=1)
    return rfp0 * rfp


def tanth_of_r(rfp_times_rfp0, opt_coeffs):
    zi = np.arctan(rfp_times_rfp0 ** 0.5 / opt_coeffs[0])
    # Minimize...
    # This is a way to minimize that
    # 1) let's us minimize the whole vector and
    # 2) doesn't use scipy, so we could do something similar in other scripting languates
    # But it's not even close to optimal.
    # If you pass in a lot of points this might fill up your memory though...
    # ----
    bounds = (0, 0.55)
    resolution = 0.001
    x = np.arange(*bounds, resolution)[:, None] * np.ones((1, len(zi)))
    for count in range(3):
        minimize = np.asarray(
            [
                (r_of_tanth(zi + x[idx], opt_coeffs) - rfp_times_rfp0) ** 2
                for idx in range(x.shape[0])
            ]
        )
        argmin = np.argmin(minimize, axis=0)
        xmin = np.asarray([x[am, idx] for idx, am in enumerate(argmin)])
        # Every iteration, scale down the offset to be narrower around the minimum
        x = (x - xmin) * 0.25 + xmin
    # ----
    return xmin + zi


def fp_optics(xyfp, opt_coeffs):
    xy = rev_az_asym(xyfp)
    rfp_times_rfp0 = np.sum(xy ** 2, axis=1) ** 0.5
    phirad = np.arctan2(-xy[:, 1], -xy[:, 0])
    phideg = phirad / (np.pi / 180.0) % 360
    thetarad = tanth_of_r(rfp_times_rfp0, opt_coeffs)
    thetadeg = thetarad / (np.pi / 180.0)
    return phideg, 90.0 - thetadeg


def sphereToCart(ras, decs):
    """Convert 3d spherical coordinates to cartesian"""
    rarads = (np.pi / 180.0) * ras
    decrads = (np.pi / 180.0) * decs
    sinras = np.sin(rarads)
    cosras = np.cos(rarads)
    sindecs = np.sin(decrads)
    cosdecs = np.cos(decrads)
    vec0s = cosras * cosdecs
    vec1s = sinras * cosdecs
    vec2s = sindecs
    return vec0s, vec1s, vec2s


def eulerm323(eul):
    mat1 = rotm1(2, eul[0])
    mat2 = rotm1(1, eul[1])
    mata = np.matmul(mat2, mat1)
    mat1 = rotm1(2, eul[2])
    rmat = np.matmul(mat1, mata)
    return rmat


def rotm1(ax, angle):
    mat = np.zeros((3, 3), dtype=np.double)
    n1 = ax
    n2 = np.mod((n1 + 1), 3)
    n3 = np.mod((n2 + 1), 3)
    sinang = np.sin(angle)
    cosang = np.cos(angle)
    mat[n1][n1] = 1.0
    mat[n2][n2] = cosang
    mat[n3][n3] = cosang
    mat[n2][n3] = sinang
    mat[n3][n2] = -sinang
    return mat


def cartToSphere(vec):
    norm = np.sqrt(np.sum(vec ** 2, axis=1))
    dec = np.arcsin(vec[:, 2] / norm)
    ra = np.arctan2(vec[:, 1], vec[:, 0])
    ra = np.mod(ra, 2.0 * np.pi)
    return ra, dec
