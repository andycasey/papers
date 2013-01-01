# This module contains a temporary Spectrum1D class for use until it's
# replaced by the upcoming astropy.specutils.Spectrum1D module.

import os
import json
import logging
import numpy as np
import pyfits


from scipy import interpolate, ndimage, polyfit, poly1d

__all__ = ['Spectrum', 'Spectrum1D', 'stitch', 'pyfits']


class Spectrum(object):
    """A temporary class holder for a multi-dimensional Spectrum object until it
    can be replaced with astropy.NDData."""

    def __init__(self, disp, flux, headers={}):
        """Initializes a `Spectrum` object with the given (multi-dimensional)
        dispersion and flux arrays.

        Inputs
        ----
        disp : `np.ndarray`
            Dispersion of the spectra.

        flux : `np.ndarray`
            Flux values for each dispersion point.

        headers : `dict`
            Headers.
        """

        self.disp = disp
        self.flux = flux
        self.headers = headers
        self.num_orders = self.flux.shape[1] if len(self.flux.shape) > 1 else len(self.flux)

        return None

    @classmethod
    def load_multispec(cls, filename):
        """Reads in a FITS multispec file into a `Spectrum` object.

        Inputs
        ----
        filename : str
            Multi-spec FITS filename to read.
        """

        if not os.path.exists(filename):
            raise IOError, "Filename '%s' does not exist." % (filename, )

        with pyfits.open(filename) as image:
            headers = image[0].header
            flux = image[0].data

        headers_dict = {}
        for k, v in headers.iteritems():

            try:
                str(v)
                json.dumps(v)

            except TypeError:
                continue

            if headers_dict.has_key(k):
                headers_dict[k] += v

            else:
                headers_dict[k] = v


        # Determine number of orders
        num_pixels = flux.shape[-1]
        num_orders = 1 if len(flux.shape) == 1 else flux.shape[-2]

        # Try linear dispersion
        try:
            crval = headers['crval1']
            crpix = headers['crpix1']
            cd = headers['cd1_1']
            ctype = headers['ctype1']

            if ctype.strip() == 'LINEAR':
                dispersion = np.zeros((num_orders, num_pixels), dtype=np.float)

                dispersion_base = (np.arange(num_pixels) + 1 - crpix) * cd + crval
                for i in xrange(num_orders):
                    dispersion[i, :] = dispersion_base

                dcflag = headers['dc-flag']
                if dcflag == 1:
                    dispersion = 10.0 ** dispersion

                elif dcflag != 0:
                    raise ValueError, "Dispersion is not linear or logarithmic (DC-FLAG = %s)" % (dcflag, )

                if num_orders == 1:
                    dispersion.shape = (num_pixels, )
                    flux = np.squeeze(flux)


                return cls(dispersion, flux, headers=headers_dict)
        
        except KeyError:
            pass

        # Get multi-spec headers
        try:
            wat = headers['wat2_*']
            num_wat_headers = len(wat)

        except KeyError:
            raise ValueError, "Cannot decipher header: need either WAT2_* or CRVAL keywords."

        # Concatenate headers
        wat_str = ''
        for i in xrange(num_wat_headers):
            # Apparently this is a hack to fix a problem in
            # older PyFits versions (< 3.1) where trailing blanks are stripped

            value = wat[i]
            if hasattr(value, 'value'): value = value.value
            value = value + (" " * (68 - len(value)))
            wat_str += value

        # Find all the spec#="..." strings
        spec_str = [''] * num_orders
        for i in xrange(num_orders):
            name = 'spec%i' % (i + 1, )

            p0 = wat_str.find(name)
            p1 = wat_str.find('"', p0)
            p2 = wat_str.find('"', p1 + 1)
            if p0 < 0 or p2 < 0 or p2 < 0:
                raise ValueError, "Cannot find '%s' in WAT2_* keyword" % (name, )

            spec_str[i] = wat_str[p1 + 1:p2]

        # Get wavelength calibration information
        w_params = np.zeros((num_orders, 9))
        w = np.zeros(9)
        for i in xrange(num_orders):
            w = np.asarray(spec_str[i].split(), dtype=np.float)
            w_params[i, :] = w
            if w[2] == -1:
                raise ValueError, "Spectrum %i has no wavelength calibration (type = %d)" % (i + 1, w[2], )

            elif w[6] != 0:
                raise ValueError, "Spectrum %i has a non-zero redshift (z=%f)" % (i + 1, w[6], )

        dispersion = np.zeros((num_orders, num_pixels), dtype=np.float)
        disp_fields = [None] * num_orders
        for i in xrange(num_orders):
            if w_params[i, 2] in (0, 1):
                
                # Simple linear or logarithmic spacing
                dispersion[i, :] = np.arange(num_pixels) * w_params[i, 4] + w_params[i, 3]
                if w_params[i, 2] == 1:
                    dispersion[i, :] = 10. ** dispersion[i, :]

            else:
                dispersion[:, i], disp_fields[i] = compute_non_linear_disp(num_pixels, spec_str[i])

        if num_orders == 1:
            flux = np.squeeze(flux)
            dispersion.shape = (num_pixels, )


        return cls(disp=dispersion, flux=flux, headers=headers_dict)





class Spectrum1D(object):
    """This is a temporary class holder for a Spectrum1D object until the
    astropy.specutils.Spectrum1D module has advanced sufficiently to replace it."""
    
    def __init__(self, disp, flux, headers={}):
        """Initializes a `Spectrum1D` object with the given dispersion and flux
        arrays.
        
        Parameters
        ----
        disp : `np.array`
            Dispersion of the spectrum (i.e. the wavelength points).
            
        flux : `np.array`
            Flux points for each `disp` point.
        """
        
        self.disp = disp
        self.flux = flux
        self.headers = headers
        
        return None
    
    @classmethod
    def load(cls, filename, **kwargs):
        """Load a Spectrum1D from a given filename.
        
        Parameters
        ----
        filename : str
            Path of the filename to load. Can be either simple FITS extension
            or an ASCII filename.
            
        Notes
        ----
        If you are loading from an non-standard ASCII file, you can pass
        kwargs to `np.loadtxt` through this function.
        """
        
        if not os.path.exists(filename):
            raise IOError("Filename '%s' does not exist." % (filename, ))
        
        if filename.endswith('.fits'):
            image = pyfits.open(filename, **kwargs)
            
            header = image[0].header
            
            if np.all([header.has_key(key) for key in ('CDELT1', 'NAXIS1', 'CRVAL1')]):
                disp = header['CRVAL1'] + np.arange(header['NAXIS1']) * header['CDELT1']
        
            disp -= header['LTV1'] if header.has_key('LTV1') else 0
            flux = image[0].data
            
            # Add the headers in
            headers = {}
            for row in header.items():
                key, value = row
                
                # Check the value is valid
                try:
                    str(value)
                    json.dumps(value)

                except TypeError:
                    continue

                if len(key) == 0 or len(str(value)) == 0: continue
                
                if headers.has_key(key):
                    if not isinstance(headers[key], list):
                        headers[key] = [headers[key]]
                    
                    headers[key].append(value)

                else:
                    headers[key] = value

            for key, value in headers.iteritems():
                if isinstance(value, list):
                    headers[key] = "\n".join(map(str, value))

                
                

        else:
            headers = {}
            disp, flux = np.loadtxt(filename, unpack=True, **kwargs)
            
        return cls(disp, flux, headers=headers)


    def save(self, filename, clobber=True):
        """Saves the `Spectrum1D` object to the specified filename.
        
        Parameters
        ----
        filename : str
            The filename to save the `Spectrum1D` object to.
            
        clobber : bool, optional
            Whether to overwite the `filename` if it already exists.
        
        Raises
        ----
        IOError
            If the filename exists and we are not asked to clobber it.
            
        ValueError
            If the ``Spectrum1D`` object does not have a linear dispersion map.
        """
        
        if os.path.exists(filename) and not clobber:
            raise IOError("Filename '%s' already exists and we have been asked not to clobber it." % (filename, ))
        
        if not filename.endswith('fits'):
            # ASCII
            
            data = np.hstack([self.disp.reshape(len(self.disp), 1), self.flux.reshape(len(self.disp), 1)])
            
            assert len(data.shape) == 2
            assert data.shape[1] == 2
            
            np.savetxt(filename, data)
            
        else:
            # FITS
            crpix1, crval1 = 1, self.disp.min()
            
            cdelt1 = np.mean(np.diff(self.disp))
            
            test_disp = (crval1 + np.arange(len(self.disp), dtype=self.disp.dtype) * cdelt1).astype(self.disp.dtype)
            
            if np.max(self.disp - test_disp) > 10e-2:
                raise ValueError("Spectrum dispersion map is not linear.")
            
            hdu = pyfits.PrimaryHDU(self.flux)
            
            headers = self.headers.copy()
            headers.update({
                'CRVAL1': crval1,
                'CRPIX1': crpix1,
                'CDELT1': cdelt1
            })
            
            for key, value in headers.iteritems():
                hdu.header.update(key, value)
            
            hdu.writeto(filename, clobber=clobber)
    
    def gaussian_smooth(self, kernel, **kwargs):
        
        # The requested kernel is in Angstroms, but the dispersion between each
        # pixel is likely less than an Angstrom, so we must calculate the true
        # kernel value
        
        true_kernel = kernel / np.median(np.diff(self.disp))
        smoothed_flux = ndimage.gaussian_filter1d(self.flux, true_kernel, **kwargs)
        
        return self.__class__(self.disp, smoothed_flux)
        

    def fit_continuum(self):


        return self


    
    def smooth(self, window_len=11, window='hanning'):
        """Smooth the data using a window with requested size.
        
        This method is based on the convolution of a scaled window with the signal.
        The signal is prepared by introducing reflected copies of the signal 
        (with the window size) in both ends so that transient parts are minimized
        in the begining and end part of the output signal.
        
        input:
            window_len: the dimension of the smoothing window; should be an odd integer
            window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
                flat window will produce a moving average smoothing.

        output:
            the smoothed spectra
            
        example:

        (This function from Heather Jacobson)
        TODO: the window parameter could be the window itself if an array instead of a string
        NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
        """

        if self.flux.size < window_len:
            raise ValueError, "Input vector needs to be bigger than window size."

        if window_len<3:
            return self

        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
            raise ValueError, "Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


        s = numpy.r_[self.flux[window_len-1:0:-1], self.flux, self.flux[-1:-window_len:-1]]
        
        if window == 'flat': #moving average
            w = numpy.ones(window_len, 'd')

        else:
            w = eval('numpy.' + window + '(window_len)', {'__builtins__': None})

        smoothed_flux = numpy.convolve(w/w.sum(), s,mode='valid')
        smoothed_flux = smoothed_flux[(window_len/2-1):-(window_len/2)]

        return self.__class__(self.disp, smoothed_flux, headers=self.headers)
    


    def doppler_correct(self, v=None, z=None, interpolate=True):
        """Performs a Doppler correction on the given `Spectrum1D` object by the
        amount required. Either velocities or redshifts can be provided.
        
        Parameters
        ----
        v : float
            The velocity (in km/s) to correct the `Spectrum1D` object by.
            
        z : float
            The redshift to correct the `Spectrum1D` object by.
            
        Raises
        ----
        ValueError
            When both a velocity `v` and a redshift `z` are provided for Doppler
            corrections.
        """
        
        if isinstance(v, float) and isinstance(z, float):
            raise ValueError('Both a velocity and redshift was supplied')
            
        c = 3e5
        
        if v is not None:
            new_disp = self.disp * np.sqrt((1+v/c)/(1-v/c))
            
        elif z is not None:
            new_disp = self.disp * (1 + z)
            
        else:
            raise TypeError('No velocity or redshift supplied.')
        
        rest_spectrum = self.__class__(new_disp, self.flux)
          
        if interpolate:
            return rest_spectrum.interpolate(self.disp)
            
        else:
            return rest_spectrum
        
    
    def fit_continuum(self, knot_spacing=200, sigma_clip=(1.0, 0.2), \
                      max_iterations=3, order=3, exclude=None, include=None, \
                      additional_points=None, function='spline', **kwargs):
        """Fits the continuum for a given `Spectrum1D` spectrum.
        
        Parameters
        ----
        knot_spacing : float or None, optional
            The knot spacing for the continuum spline function in Angstroms. Optional.
            If not provided then the knot spacing will be determined automatically.
        
        sigma_clip : a tuple of two floats, optional
            This is the lower and upper sigma clipping level respectively. Optional.
            
        max_iterations : int, optional
            Maximum number of spline-fitting operations.
            
        order : int, optional
            The order of the spline function to fit.
            
        exclude : list of tuple-types containing floats, optional
            A list of wavelength regions to always exclude when determining the
            continuum. Example:
            
            >> exclude = [
            >>    (3890.0, 4110.0),
            >>    (4310.0, 4340.0)
            >>  ]
            
            In the example above the regions between 3890 A and 4110 A, as well as
            4310 A to 4340 A will always be excluded when determining the continuum
            regions.

        function: only 'spline' or 'poly'
            
        include : list of tuple-types containing floats, optional
            A list of wavelength regions to always include when determining the
            continuum.
        """
        
        exclusions = []
        continuum_indices = range(len(self.flux))
        
        # See if there are any regions we need to exclude
        if exclude is not None and len(exclude) > 0:
            exclude_indices = []
            
            if isinstance(exclude[0], float) and len(exclude) == 2:
                # Only two floats given, so we only have one region to exclude
                exclude_indices.extend(range(*np.searchsorted(self.disp, exclude)))
                
            else:
                # Multiple regions provided
                for exclude_region in exclude:
                    exclude_indices.extend(range(*np.searchsorted(self.disp, exclude_region)))
        
            continuum_indices = np.sort(list(set(continuum_indices).difference(exclude_indices)))
        
        # See if there are any regions we should always include
        if include is not None and len(include) > 0:
            include_indices = []
            
            if isinstance(include[0], float) and len(include) == 2:
                # Only two floats given, so we can only have one region to include
                include_indices.extend(range(*np.searchsorted(self.disp, include)))
                
            else:
                # Multiple regions provided
                for include_region in include:
                    include_indices.extend(range(*np.searchsorted(self.disp, include_region)))
        

        # We should exclude non-finite numbers from the fit
        non_finite_indices = np.where(~np.isfinite(self.flux))[0]
        continuum_indices = np.sort(list(set(continuum_indices).difference(non_finite_indices)))

        # We should also exclude zero or negative flux points from the fit
        zero_flux_indices = np.where(0 >= self.flux)[0]
        continuum_indices = np.sort(list(set(continuum_indices).difference(zero_flux_indices)))

        if knot_spacing is None or knot_spacing == 0:
            knots = []

        else:
            knot_spacing = abs(knot_spacing)
            
            end_spacing = ((self.disp[-1] - self.disp[0]) % knot_spacing) /2.
        
            if knot_spacing/2. > end_spacing: end_spacing += knot_spacing/2.
                
            knots = np.arange(self.disp[0] + end_spacing, self.disp[-1] - end_spacing + knot_spacing, knot_spacing)
            if len(knots) > 0 and knots[-1] > self.disp[continuum_indices][-1]:
                knots = knots[:knots.searchsorted(self.disp[continuum_indices][-1])]
                
            if len(knots) > 0 and knots[0] < self.disp[continuum_indices][0]:
                knots = knots[knots.searchsorted(self.disp[continuum_indices][0]):]


        for iteration in xrange(max_iterations):
            
            splrep_disp = self.disp[continuum_indices]
            splrep_flux = self.flux[continuum_indices]

            splrep_weights = np.ones(len(splrep_disp))

            # We need to add in additional points at the last minute here
            if additional_points is not None and len(additional_points) > 0:

                for point, flux, weight in additional_points:

                    # Get the index of the fit
                    insert_index = np.searchsorted(splrep_disp, point)
                    
                    # Insert the values
                    splrep_disp = np.insert(splrep_disp, insert_index, point)
                    splrep_flux = np.insert(splrep_flux, insert_index, flux)
                    splrep_weights = np.insert(splrep_weights, insert_index, weight)



            if function == 'spline':
                order = 5 if order > 5 else order
                

                tck = interpolate.splrep(splrep_disp, \
                                         splrep_flux, \
                                         k=order, task=-1, t=knots, w=splrep_weights)

                continuum = interpolate.splev(self.disp, tck)

            elif function == 'poly':
            
                p = poly1d(polyfit(splrep_disp, splrep_flux, order))
                continuum = p(self.disp)

            else:
                raise ValueError("Unknown function type: only spline or poly available")
            

            difference = continuum - self.flux
        
            #n = 100./np.median(np.diff(self.disp))
            #n = 100
            #convolved_flux = np.convolve(self.flux, np.ones(n)/n)
            #sigma_difference = (difference - np.median(difference)) / np.std(difference)
            sigma_difference = difference / np.std(difference)

            # Clip 
            upper_exclude = np.where(sigma_difference > sigma_clip[1])[0]
            lower_exclude = np.where(sigma_difference < -sigma_clip[0])[0]
            
            exclude_indices = list(upper_exclude)
            exclude_indices.extend(lower_exclude)
            exclude_indices = np.array(exclude_indices)
            
            if len(exclude_indices) is 0: break
            
            exclusions.extend(exclude_indices)
            
            # Before excluding anything, we must check to see if there are regions
            # which we should never exclude
            if include is not None:
                exclude_indices = set(exclude_indices).difference(include_indices)
            
            # Remove regions that have been excluded
            continuum_indices = np.sort(list(set(continuum_indices).difference(exclude_indices)))
        
        
        return self.__class__(disp=self.disp, flux=continuum, headers=self.headers)
        


    
    def interpolate(self, new_disp, mode='linear', bounds_error=False,
                    fill_value=np.nan):
        """Interpolate the `Spectrum1D` onto a new dispersion map.
        
        Parameters
        ----
        new_disp : np.array
            An array of floating-point types containing the new dispersion points.
            
        mode : str
            Interpolation mode. See `scipy.interpolate.interp1d` for available
            options.
        
        bounds_error : bool
            See `scipy.interpolate.interp1d` for details.
        
        fill_value : float-type
            See `scipy.interpolate.interp1d`
        """
        
        f = interpolate.interp1d(self.disp, self.flux, kind=mode, copy=False,
                                 bounds_error=bounds_error, fill_value=fill_value)
        
        if isinstance(new_disp, float):
            return self.__class__(new_disp, f(new_disp))
            
        else:
            return self.__class__(new_disp, f(new_disp))
        
        
        

def fxcor(spectrum, template, wl_start, wl_end):
    
    try:
        # [TODO] Destroy IRAF
        from pyraf import iraf
        
    except NameError:
        iraf.noao()
        iraf.noao.onedspec()
        iraf.rv()
    
    if not os.path.exists(spectrum):
        raise IOError('Spectrum file "%s" does not exist.' % (spectrum, ))
    
    if not os.path.exists(template):
        raise IOError('Template file "%s" does not exist.' % (template, ))
        
    
    # Check whether the spectrum and template are ASCII files. If they are we
    # need them into FITS to continue.
    to_remove = [] # for cleaning up later
    
    spectrum_base, spectrum_ext = os.path.splitext(spectrum)
    spectrum_ext = spectrum_ext[1:]
    if spectrum_ext not in ('txt', 'ascii', 'fits'):
        logging.warn("Don't know what format the spectrum provided is in. Assuming ASCII..")
    
    if spectrum_ext != 'fits':
        wave = np.loadtxt(spectrum, usecols=(0, ))
        
        spectrum_fits = '.'.join([spectrum_base, '_tmp.fits'])
        iraf.onedspec.rspectext(input = spectrum,
                                output = spectrum_fits,
                                dtype = 'linear',
                                crval1 = np.min(wave),
                                cdelt1 = np.median(np.diff(wave)))
        to_remove.append(spectrum_fits)
    else: spectrum_fits = spectrum
    
    template_base, template_ext = os.path.splitext(template)
    template_ext = template_ext[1:]
    if template_ext not in ('txt', 'ascii', 'fits'):
        logging.warn("Don't know what format the template provided is in. Assuming ASCII..")
    
    if template_ext != 'fits':
        wave = np.loadtxt(spectrum, usecols=(0, ))
        template_fits = '.'.join([template_base, '_tmp.fits'])
        iraf.onedspec.rspectext(input = template,
                                output = template_fits,
                                dtype = 'linear',
                                crval1 = np.min(wave),
                                cdelt1 = np.median(np.diff(wave)))
        to_remove.append(template_fits)
    else: template_fits = template
    
    sample = ':'.join(map(str, [wl_start, wl_end]))
    output_filename = 'fxcor'
    
    try:
    
        iraf.rv.fxcor(objects=spectrum_fits, templates=template_fits, rebin='template',
            osample=sample, rsample=sample, output=output_filename, verbose='long', interactive=0, Stdout=1)
    
    except:
        raise
    
    else:
        
        fxcor_log = []
        with open(output_filename + '.log', 'r') as output:
            fxcor_log.extend(output.readlines())
        
        fxcor_out = []
        with open(output_filename + '.txt', 'r') as output:
            fxcor_out.extend(output.readlines())
        
        fxcor_log = ''.join(fxcor_log)
        fxcor_out = ''.join(fxcor_out)
        
        fxcor_output = np.loadtxt(output_filename + '.txt', comments='#', dtype=str)
        fxcor_hght, fxcor_fwhm, fxcor_tdr, fxcor_obs, fxcor_rel, fxcor_helio, fxcor_err = fxcor_output[-7::]
            
        # Translate
        if 'INDEF' in [fxcor_hght, fxcor_rel, fxcor_err]:
            v_rad = fxcor_rel.replace('INDEF', '0')
            correlation = fxcor_hght.replace('INDEF', '0')
            v_err = fxcor_rel.replace('INDEF', '-999')
            
        v_rad, correlation, v_err = map(float, [fxcor_rel, fxcor_hght, fxcor_err])
    
        os.system('rm -f %s fxcor.txt fxcor.log fxcor.gki' % (' '.join(to_remove), ))
    
        return (v_rad, v_err, correlation, )




def compute_non_linear_disp(nwave, specstr, verbose=False):
    """Compute non-linear wavelengths from multispec string
    
    Returns wavelength array and dispersion fields.
    Raises a ValueError if it can't understand the dispersion string.
    """

    fields = specstr.split()
    if int(fields[2]) != 2:
        raise ValueError('Not nonlinear dispersion: dtype=' + fields[2])
    if len(fields) < 12:
        raise ValueError('Bad spectrum format (only %d fields)' % len(fields))
    wt = float(fields[9])
    w0 = float(fields[10])
    ftype = int(fields[11])
    if ftype == 3:

        # cubic spline

        if len(fields) < 15:
            raise ValueError('Bad spline format (only %d fields)' % len(fields))
        npieces = int(fields[12])
        pmin = float(fields[13])
        pmax = float(fields[14])
        if verbose:
            print 'Dispersion is order-%d cubic spline' % npieces
        if len(fields) != 15+npieces+3:
            raise ValueError('Bad order-%d spline format (%d fields)' % (npieces,len(fields)))
        coeff = np.asarray(fields[15:],dtype=float)
        # normalized x coordinates
        s = (np.arange(nwave,dtype=float)+1-pmin)/(pmax-pmin)*npieces
        j = s.astype(int).clip(0, npieces-1)
        a = (j+1)-s
        b = s-j
        x0 = a**3
        x1 = 1+3*a*(1+a*b)
        x2 = 1+3*b*(1+a*b)
        x3 = b**3
        wave = coeff[j]*x0 + coeff[j+1]*x1 + coeff[j+2]*x2 + coeff[j+3]*x3

    elif ftype == 1 or ftype == 2:

        # chebyshev or legendre polynomial
        # legendre not tested yet

        if len(fields) < 15:
            raise ValueError('Bad polynomial format (only %d fields)' % len(fields))
        order = int(fields[12])
        pmin = float(fields[13])
        pmax = float(fields[14])
        if verbose:
            if ftype == 1:
                print 'Dispersion is order-%d Chebyshev polynomial' % order
            else:
                print 'Dispersion is order-%d Legendre polynomial (NEEDS TEST)' % order
        if len(fields) != 15+order:
            raise ValueError('Bad order-%d polynomial format (%d fields)' % (order, len(fields)))
        coeff = np.asarray(fields[15:],dtype=float)
        # normalized x coordinates
        pmiddle = (pmax+pmin)/2
        prange = pmax-pmin
        x = (np.arange(nwave,dtype=float)+1-pmiddle)/(prange/2)
        p0 = np.ones(nwave,dtype=float)
        p1 = x
        wave = p0*coeff[0] + p1*coeff[1]
        for i in range(2, order):
            if ftype == 1:
                # chebyshev
                p2 = 2*x*p1 - p0
            else:
                # legendre
                p2 = ((2*i-1)*x*p1-(i-1)*p0) / i
            wave = wave + p2*coeff[i]
            p0 = p1
            p1 = p2

    else:
        raise ValueError('Cannot handle dispersion function of type %d' % ftype)

    return wave, fields


def stitch(spectra, wl_ranges=None, mode='average'):
    """Stitches together two overlapping spectra.

    Inputs
    ----
    spectra : list of `Spectrum1D` objects.
        A list of spectra to stitch.

    wl_ranges : list of tuples with length two, optional
        The wavelength regions to use for stitching each spectrum. If no wavelength
        region is specified then the entire spectrum is used for stitching. If only
        wavelength region for some of the spectra are given, then None should be
        supplied when to use the entire spectrum.

    mode : str
        How to stack the spectra on a pixel-by-pixel basis. Available: average, sum, quadrature

    Returns
    ----
    stitched_spectrum : `Spectrum1D`
        A linear stitched spectrum.
    """

    if mode not in ('average', 'sum', 'quadrature', ):
        raise ValueError, "Mode must be either 'average', 'quadrature' or 'sum'."


    # Ensure that our spectra go from blue to red
    wl_starts = []
    wl_ends = []
    min_disp_step = 999
    for spectrum in spectra:
        wl_starts.append(spectrum.disp[0])
        wl_ends.append(spectrum.disp[-1])

        if np.min(np.diff(spectrum.disp)) < min_disp_step:
            min_disp_step = np.min(np.diff(spectrum.disp))

    sorted_indices = np.argsort(wl_starts)

    # Let's generate our new dispersion map
    new_disp = np.arange(spectra[sorted_indices[0]].disp[0], spectra[sorted_indices[-1]].disp[-1], min_disp_step)
    new_flux = np.zeros(len(new_disp))
    
    # And an array with the number of contributing flux points for each pixel point
    num_flux = np.zeros(len(new_disp))

    # Maintain headers, give preference to blue
    headers = {}
    for index in sorted_indices[::-1]:
        headers.update(spectra[index].headers)
        

    for index in sorted_indices:
        spectrum = spectra[index]


        spectrum.flux[np.where(~np.isfinite(spectrum.flux))] = 0.0

        if wl_ranges is None or wl_ranges[index] is None:
            wl_range = [spectrum.disp[0], spectrum.disp[-1]]

        else:
            wl_range = wl_ranges[index]
            

        # Get the subset dispersion map
        dispersion_indices = np.searchsorted(new_disp, wl_range)
        subset_disp = new_disp[dispersion_indices[0]:dispersion_indices[1]]

        # Interpolate the spectra
        f = interpolate.interp1d(spectrum.disp, spectrum.flux, kind='linear', copy=False,
                                 bounds_error=False, fill_value=np.nan)

        subset_flux = f(subset_disp)

        # Put the subset flux into our master arrays
        additional_flux = pow(subset_flux, 2) if mode == 'quadrature' else subset_flux

        good_indices = np.where(additional_flux > 0.)[0]

        new_flux[good_indices + dispersion_indices[0]] += additional_flux[good_indices]
        num_flux[good_indices + dispersion_indices[0]] += 1


        
    
    if mode == 'average':
        new_flux /= num_flux

    elif mode == 'quadrature':
        new_flux = pow(new_flux, 0.5)

    # Remove out non-finite values on the edge
    idx_l = np.searchsorted(np.isfinite(new_flux), True)
    idx_r = len(new_flux) - np.searchsorted(np.isfinite(new_flux)[::-1], True)

    return Spectrum1D(new_disp[idx_l:idx_r], new_flux[idx_l:idx_r], headers=headers)

