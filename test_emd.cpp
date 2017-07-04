
#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <list>
#include <string>

#include "H5Cpp.h"
#include "hdf5.h"

using namespace H5;

#define MAX_NAME 1024

void do_dtype(hid_t);
void do_dset(hid_t);
void do_link(hid_t, char *);
void scan_group(hid_t);
void do_attr(hid_t);
void scan_attrs(hid_t);
void do_plist(hid_t);

class EMDObject {
  protected:
    std::string name;
    int i;
    ssize_t len;
    hsize_t nobj;
    herr_t err;
    int otype;
    hid_t type_id;
    hid_t dsid;

    // datatype
    H5T_class_t t_class;
    // size of datatype
    size_t t_size;
    // Integer data can be signed two's complement (H5T_SGN_2) or unsigned (H5T_SGN_NONE)
    H5T_sign_t int_sign;

  public:
    /*
     *  Analyze a data type description
     */
    void
      do_dtype(hid_t tid) {
        t_class = H5Tget_class(tid);
        t_size = H5Tget_size(tid);
        printf(" Datatype size %d.\n", t_size);

        if(t_class < 0){
          //puts(" Invalid datatype.\n");
        } else {
          /*
           * Each class has specific properties that can be
           * retrieved, e.g., size, byte order, exponent, etc.
           */
          if(t_class == H5T_INTEGER) {
            puts(" Datatype is 'H5T_INTEGER'.\n");
            int_sign = H5Tget_sign(tid);
            std::string _unsigned_text = (int_sign == H5T_SGN_NONE) ? " unsigned int " : " signed int ";
            std::cout << "##### SIGNAL: " <<  _unsigned_text << std::endl;
            std::cout << " sizeof(int): " << sizeof(int) << std::endl;
            std::cout << " sizeof(unsigned short): " << sizeof(unsigned short) << std::endl;

             /* display size, signed, endianess, etc. */
          } else if(t_class == H5T_FLOAT) {
            //puts(" Datatype is 'H5T_FLOAT'.\n");
            /* display size, endianess, exponennt, etc. */
          } else if(t_class == H5T_STRING) {
            //puts(" Datatype is 'H5T_STRING'.\n");
            /* display size, padding, termination, etc. */
          } else if(t_class == H5T_BITFIELD) {
            //puts(" Datatype is 'H5T_BITFIELD'.\n");
            /* display size, label, etc. */
          } else if(t_class == H5T_OPAQUE) {
            //puts(" Datatype is 'H5T_OPAQUE'.\n");
            /* display size, etc. */
          } else if(t_class == H5T_COMPOUND) {
            //puts(" Datatype is 'H5T_COMPOUND'.\n");
            /* recursively display each member: field name, type  */
          } else if(t_class == H5T_ARRAY) {
            //puts(" Datatype is 'H5T_COMPOUND'.\n");
            /* display  dimensions, base type  */
          } else if(t_class == H5T_ENUM) {
            //puts(" Datatype is 'H5T_ENUM'.\n");
            /* display elements: name, value   */
          } else  {
            //puts(" Datatype is 'Other'.\n");
            /* eg. Object Reference, ...and so on ... */
          }
        }
      }

    /*
     *  Analyze a symbolic link
     *
     * The main thing you can do with a link is find out
     * what it points to.
     */
    void
      do_link(hid_t gid, char *name) {
        herr_t status;
        char target[MAX_NAME];

        status = H5Gget_linkval(gid, name, MAX_NAME, target  ) ;
        printf("Symlink: %s points to: %s\n", name, target);
      }
};

class EMDAttribute : public EMDObject {
hid_t aid;
ssize_t len;
hid_t atype;
hid_t aspace;
public:
  EMDAttribute( hid_t id ){
    aid = id;
  }

  /*
   *  Process one attribute.
   *  This is similar to the information about a dataset.
   */
  void do_attr(hid_t aid) {

    char buf[MAX_NAME];
    /*
     * Get the name of the attribute.
     */
    len = H5Aget_name(aid, MAX_NAME, buf );
    printf("Attribute Name : %s\n",buf);
    name = buf;

    /*
     * Get attribute information: dataspace, data type
     */
    aspace = H5Aget_space(aid); /* the dimensions of the attribute data */
    atype  = H5Aget_type(aid);
    do_dtype(atype);

    /*
     * The datatype and dataspace can be used to read all or
     * part of the data.  (Not shown in this example.)
     */

    /* ... read data with H5Aread, write with H5Awrite, etc. */

    H5Tclose(atype);
    H5Sclose(aspace);
  }

};

class EMDDataSet : public EMDObject {

private:
  std::vector<EMDAttribute*> _attributes;

  // size of chunks for the raw data of a chunked layout dataset
  hsize_t chunk_dims_out[3];
  int  rank_chunk;
  int nfilters;

  char* raw_data;
  size_t size_raw_data;

  public:
    EMDDataSet(   hid_t id ){
      dsid = id;
    }

    /*
     *  Run through all the attributes of a dataset or group.
     *  This is similar to iterating through a group.
     */

    void scan_attrs(hid_t oid) {
        int na;
        hid_t aid;
        int i;

        na = H5Aget_num_attrs(oid);

        for (i = 0; i < na; i++) {
          aid =	H5Aopen_idx(oid, (unsigned int)i );
          EMDAttribute* new_attribute = new EMDAttribute( aid );
          new_attribute->do_attr( aid );
          _attributes.push_back( new_attribute );
          H5Aclose(aid);
        }
      }

    /*
     *  Retrieve information about a dataset.
     *
     *  Many other possible actions.
     *
     *  This example does not read the data of the dataset.
     */
    void do_dset(hid_t did){
        hid_t tid;
        hid_t pid;
        // dataset space
        hid_t space_id;
        hsize_t size;
        char ds_name[MAX_NAME];


        /*
         * Information about the group:
         *  Name and attributes
         *
         *  Other info., not shown here: number of links, object id
         */
        H5Iget_name(did, ds_name, MAX_NAME  );
        name = ds_name;

        //printf("Dataset Name : ");
        //puts(ds_name);
        //printf("\n");
        //printf(" ATTRIBUTES:\n");
        /*
         *  process the attributes of the dataset, if any.
         */
        scan_attrs(did);
        /*
         * Get dataset information: dataspace, data type
         */
        space_id = H5Dget_space(did); /* the dimensions of the dataset (not shown) */
        tid = H5Dget_type(did);
        //printf(" DATA TYPE:\n");
        do_dtype(tid);

        /*
         * Retrieve and analyse the dataset properties
         */
        //printf(" PROPERTY LIST:\n");
        pid = H5Dget_create_plist(did); /* get creation property list */
        do_plist(pid);
        //printf(" STORAGE LIST:\n");
        size = H5Dget_storage_size(did);
        printf("Total space currently written in file: %d\n\n\n",(int)size);

        size_raw_data = chunk_dims_out[0];
        for(size_t pos_rank = 1; pos_rank < rank_chunk; pos_rank++ ){
          size_raw_data *= chunk_dims_out[pos_rank];
        }

        unsigned short * dset_data = new unsigned short[size_raw_data];
        size_raw_data *= t_size;
        std::cout << " size_raw_data " << size_raw_data << std::endl;

        herr_t      status;

        /*
         * Create a datatype to refer to.
         */
         hid_t		type_id;       /* Datatype ID	   	 */
         hid_t		size1;       /* Datatype ID	   	 */

      type_id = H5Tvlen_create (H5T_NATIVE_USHORT);
      size1 = H5Dget_storage_size(did);
      std::vector<char> buf(static_cast<int>(size1), 0x00); //Allocate and Zero the array
      type_id = H5Dget_type(did);
      status = H5Dread(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(buf.front()) );
      std::cout << " buf.size " << buf.size() << std::endl;



        /*
         * The datatype and dataspace can be used to read all or
         * part of the data.  (Not shown in this example.)
         */

        /* ... read data with H5Dread, write with H5Dwrite, etc. */
        H5Pclose(pid);
        H5Tclose(tid);
        H5Sclose(space_id);
      }

      /*
       *   Example of information that can be read from a Dataset Creation
       *   Property List.
       *
       *   There are many other possibilities, and there are other property
       *   lists.
       */
      void do_plist(hid_t pid) {

          H5Z_filter_t  filtn;
          int i;
          unsigned int   filt_flags, filt_conf;
          size_t cd_nelmts;
          unsigned int cd_values[32] ;
          char f_name[MAX_NAME];
          H5D_fill_time_t ft;
          H5D_alloc_time_t at;
          H5D_fill_value_t fvstatus;
          unsigned int szip_options_mask;
          unsigned int szip_pixels_per_block;

          /* zillions of things might be on the plist */
          /*  here are a few... */

          /*
           * get chunking information: rank and dimensions.
           *
           *  For other layouts, would get the relevant information.
           */
          if(H5D_CHUNKED == H5Pget_layout(pid)){
            rank_chunk = H5Pget_chunk(pid, 3, chunk_dims_out);
            printf("chunk rank %d, dimensions %lu x %lu x %lu\n", rank_chunk,
                (unsigned long)(chunk_dims_out[0]),
                (unsigned long)(chunk_dims_out[1]),
              (unsigned long)(chunk_dims_out[2]));
          } /* else if contiguous, etc. */

          /*
           *  Get optional filters, if any.
           *
           *  This include optional checksum and compression methods.
           */

          nfilters = H5Pget_nfilters(pid);
          for (i = 0; i < nfilters; i++)
          {
            /* For each filter, get
             *   filter ID
             *   filter specific parameters
             */
            cd_nelmts = 32;
            filtn = H5Pget_filter(pid, (unsigned)i,
                &filt_flags, &cd_nelmts, cd_values,
                (size_t)MAX_NAME, f_name, &filt_conf);
            /*
             *  These are the predefined filters
             */
            switch (filtn) {
              case H5Z_FILTER_DEFLATE:  /* AKA GZIP compression */
                printf("DEFLATE level = %d\n", cd_values[0]);
                break;
              case H5Z_FILTER_SHUFFLE:
                printf("SHUFFLE\n"); /* no parms */
                break;
              case H5Z_FILTER_FLETCHER32:
                printf("FLETCHER32\n"); /* Error Detection Code */
                break;
              case H5Z_FILTER_SZIP:
                szip_options_mask=cd_values[0];;
                szip_pixels_per_block=cd_values[1];

                printf("SZIP COMPRESSION: ");
                printf("PIXELS_PER_BLOCK %d\n",
                    szip_pixels_per_block);
                /* print SZIP options mask, etc. */
                break;
              default:
                printf("UNKNOWN_FILTER\n" );
                break;
            }
          }

          /*
           *  Get the fill value information:
           *    - when to allocate space on disk
           *    - when to fill on disk
           *    - value to fill, if any
           */

          printf("ALLOC_TIME ");
          H5Pget_alloc_time(pid, &at);

          switch (at)
          {
            case H5D_ALLOC_TIME_EARLY:
              printf("EARLY\n");
              break;
            case H5D_ALLOC_TIME_INCR:
              printf("INCR\n");
              break;
            case H5D_ALLOC_TIME_LATE:
              printf("LATE\n");
              break;
            default:
              printf("unknown allocation policy");
              break;
          }

          printf("FILL_TIME: ");
          H5Pget_fill_time(pid, &ft);
          switch ( ft )
          {
            case H5D_FILL_TIME_ALLOC:
              printf("ALLOC\n");
              break;
            case H5D_FILL_TIME_NEVER:
              printf("NEVER\n");
              break;
            case H5D_FILL_TIME_IFSET:
              printf("IFSET\n");
              break;
            default:
              printf("?\n");
              break;
          }

          H5Pfill_value_defined(pid, &fvstatus);

          if (fvstatus == H5D_FILL_VALUE_UNDEFINED)
          {
            printf("No fill value defined, will use default\n");
          } else {
            /* Read  the fill value with H5Pget_fill_value.
             * Fill value is the same data type as the dataset.
             * (details not shown)
             **/
          }

          /* ... and so on for other dataset properties ... */
        }
};

class EMDGroup : public EMDObject {
  private:
    hid_t gid;
    hid_t grpid;
    char group_name[MAX_NAME];
    char memb_name[MAX_NAME];
    std::vector<EMDGroup*> _groups;
    std::vector<EMDDataSet*> _datasets;
    std::vector<EMDAttribute*> _attributes;

  public:
    EMDGroup( hid_t id ){
      gid = id;
    }

    /*
     *  Run through all the attributes of a dataset or group.
     *  This is similar to iterating through a group.
     */
    void scan_attrs(hid_t oid) {
        int na;
        hid_t aid;
        int i;

        na = H5Aget_num_attrs(oid);

        for (i = 0; i < na; i++) {
          aid =	H5Aopen_idx(oid, (unsigned int)i );
          EMDAttribute* new_attribute = new EMDAttribute( aid );
          new_attribute->do_attr( aid );
          _attributes.push_back( new_attribute );
          H5Aclose(aid);
        }
      }


    /*
     * Process a group and all it's members
     *
     *   This can be used as a model to implement different actions and
     *   searches.
     */

    void scan_group(hid_t gid) {

      /*
       * Information about the group:
       *  Name and attributes
       *
       *  Other info., not shown here: number of links, object id
       */
      len = H5Iget_name (gid, group_name, MAX_NAME);
      //printf("Group Name: %s\n",group_name);
      name = group_name;

      /*
       *  process the attributes of the group, if any.
       */
      scan_attrs(gid);

      /*
       *  Get all the members of the groups, one at a time.
       */
      err = H5Gget_num_objs(gid, &nobj);
      for (i = 0; i < nobj; i++) {
        /*
         *  For each object in the group, get the name and
         *   what type of object it is.
         */
        //printf("  Member #: %d ",i);fflush(stdout);
        len = H5Gget_objname_by_idx(gid, (hsize_t)i, memb_name, (size_t)MAX_NAME );
        //printf("  length %d ",len);
        //fflush(stdout);
        //printf("  Member: %s ",memb_name);
        //fflush(stdout);
        otype = H5Gget_objtype_by_idx(gid, (size_t)i );
        /*
         * process each object according to its type
         */
        switch(otype) {
          case H5G_LINK:
            {
              //printf(" SYM_LINK:\n");
              do_link(gid,memb_name);
              break;
            }
          case H5G_GROUP:
            {
              //printf(" GROUP:\n");
              grpid = H5Gopen(gid,memb_name, H5P_DEFAULT);
              EMDGroup* new_group = new EMDGroup( grpid );
              new_group->scan_group( grpid );
              _groups.push_back( new_group );
              H5Gclose(grpid);
              break;
            }
          case H5G_DATASET:
            {
              //printf(" DATASET:\n");
              dsid = H5Dopen(gid,memb_name, H5P_DEFAULT);
              EMDDataSet* new_dataset = new EMDDataSet( dsid );
              new_dataset->do_dset( dsid );
              _datasets.push_back( new_dataset );
              //do_dset(dsid);
              H5Dclose(dsid);
              break;
            }
          case H5G_TYPE:
            {
              //printf(" DATA TYPE:\n");
              type_id = H5Topen(gid,memb_name, H5P_DEFAULT);
              do_dtype(type_id);
              H5Tclose(type_id);
              break;
            }
          default:
            {
              //printf(" unknown?\n");
              break;
            }
        }
      }
    }

};

const H5std_string FILE_NAME( "../STEM_2016-12-07_14h43m22s.emd" );
const H5std_string DATASET_NAME( "Data" );

int main(int argc, char **argv){
  bool valid_file = H5::H5File::isHdf5( FILE_NAME );
  std::cout << " hdf5 format? " << std::boolalpha << valid_file << std::endl;
  if (valid_file ){
    try {
      hid_t    grp_id;
      herr_t   status;
      /*
       *  Example: open a file, open the root, scan the whole file.
       */
      H5File *file = new H5File( FILE_NAME, H5F_ACC_RDWR );
      grp_id = H5Gopen(file->getId(),"/", H5P_DEFAULT);
      EMDGroup grp(grp_id);
      grp.scan_group( grp_id );
      status = H5Fclose(file->getId());
    } catch(const H5::FileIException&) {
    }
  }
  return 0;
}
