use nom::{bytes::complete::take, IResult, number::complete::{be_u16, be_u32}};

// Each snapgene packet has a header (5 bytes, tag + length of packet) and then the data which can be of a few packet "types"
#[derive(Debug, PartialEq)]
struct SnapgenePacket {
    packet_header: PacketHeader,
    packet_data: PacketTypes,
}

#[derive(Debug, PartialEq)]
enum PacketTypes {
    Cookie(CookiePacket),
    Sequence(SequencePacket),
    // Notes(NotesPacket),
    // Feature(FeaturePacket),
    // Primer(PrimerPacket),
}

// A packet header is the first 5 bytes of a packet
#[derive(Debug, PartialEq)]
struct PacketHeader {
    packet_tag: PacketType,
    packet_length: u32,
}

// Cookie is the initial packet
#[derive(Debug, PartialEq)]
struct CookiePacket {
    packet_header: PacketHeader,
    magic_cookie: [u8; 8],
    seq_type: SeqType,
    export_version: u16,
    import_version: u16,
}

// Appears once
#[derive(Debug, PartialEq)]
struct SequencePacket {
    packet_header: PacketHeader,
    sequence_flag: Topology,
    seq: Vec<u8>,
}

#[derive(Debug, PartialEq)]
enum Topology {
    Linear = 0x00, // Absent flag = linear
    Circular,      // Otherwise circular
}

#[derive(Debug, PartialEq)]
enum SeqType {
    DNA = 1,
    RNA = 2,
    Protein = 3,
}

#[derive(Debug, PartialEq)]
enum PacketType {
    Cookie = 0x09,
    Sequence = 0x00,
    Notes = 0x06,
    Feature = 0x0A,
}


impl PacketHeader {
    fn parse_header(input: &[u8]) -> IResult<&[u8], PacketHeader> {
        // Header is tag (1 byte) + length (4 bytes) in big endian (long)
        let (input, packet_tag) = take(1usize)(input)?;
        let (input, packet_length) = be_u32(input)?;

        let packet_header = PacketHeader {
            packet_tag: match packet_tag[0] {
                0x00 => PacketType::Sequence,
                0x06 => PacketType::Notes,
                0x09 => PacketType::Cookie,
                0x0A => PacketType::Feature,
                _ => panic!("Unknown packet type: byte = {:?}", packet_tag[0]),
            },
            // Convert to u32
            packet_length: packet_length,
        };

        Ok((input, packet_header))
    }
}


// This is a binary parser
impl CookiePacket {
    fn parse_cookie(input: &[u8]) -> IResult<&[u8], CookiePacket> {
        // Cookie always starts with 0x09, return otherwise
        if input[0] != 0x09 {
            panic!("Incorrect file or corrupted snapgene file");
        }

        // Parse the packet header (first 5 bytes), and then the rest of the cookie
        let (input, packet_header) = PacketHeader::parse_header(input)?;
        let (input, magic_cookie) = take(8usize)(input)?;
        let (input, seq_type) = take(2usize)(input)?;
        // let (input, export_version) = take(2usize)(input)?;
        // let (input, import_version) = take(2usize)(input)?;
        let (input, export_version) = be_u16(input)?;
        let (input, import_version) = be_u16(input)?;

        // Set the values
        let cookie_packet = CookiePacket {
            packet_header,
            magic_cookie: magic_cookie.try_into().unwrap(),
            seq_type: match u16::from_be_bytes(seq_type.try_into().unwrap()) {
                1 => SeqType::DNA,
                2 => SeqType::RNA,
                3 => SeqType::Protein,
                _ => panic!("Unknown sequence type"),
            },
            export_version,
            import_version
        };

        Ok((input, cookie_packet))
    }
}

impl SequencePacket {
    fn parse_dnapacket(input: &[u8]) -> IResult<&[u8], SequencePacket> {
        // Parse the packet header
        let (input, packet_header) = PacketHeader::parse_header(input)?;
        let (input, sequence_flag) = take(1usize)(input)?;
        // Sequence length is packet length - 1 (due to sequence flag)
        let length = packet_header.packet_length - 1;
        let (input, seq) = take(length as usize)(input)?;

        // Set the values
        let snap_parser = SequencePacket {
            packet_header,
            sequence_flag: match sequence_flag[0] {
                0x00 => Topology::Linear,
                _ => Topology::Circular,
            },
            seq: seq.to_vec(),
        };

        Ok((input, snap_parser))
    }
}


fn main() {
    let f = std::fs::read("./pcDNA4 TO.dna").unwrap();
    // let (input, cookie_packet) = CookiePacket::parse_cookie(&f).unwrap();
    // let (input, sequence) = SequencePacket::parse_dnapacket(input).unwrap();
    // let sequence = String::from_utf8(sequence.seq).unwrap();

    // Repeatedly parse the input until we get to the end, using the various parsers


}
